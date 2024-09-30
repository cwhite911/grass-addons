#!/usr/bin/env python3
#
############################################################################
#
# MODULE:	r.slopeunits for GRASS 7
# AUTHOR(S):    Ivan Marchesini, Massimiliano Alvioli
# PURPOSE:	To create a raster layer of slope units
# COPYRIGHT: (C) 2004-2012 by the GRASS Development Team
#
# 		This program is free software under the GNU General Public
# 		License (>=v2). Read the file COPYING that comes with GRASS
# 		for details.
#
#############################################################################
#
# %module
# % description: Create a raster layer of slope units
# % keywords: raster
# % keywords: elevation
# % keywords: slopeunits
# %end

# %option G_OPT_R_INPUT
# % key: demmap
# % description: Input digital elevation model
# % required : yes
# %end

# %option G_OPT_R_INPUT
# % key: plainsmap
# % description: Input raster map of alluvial_plains
# % required : no
# %end

# %option G_OPT_R_OUTPUT
# % key: slumap
# % description: Output Slope Units layer (the main output)
# % required : yes
# %end

# %option G_OPT_R_OUTPUT
# % key: slumapclean
# % description: Output Slope Units layer (the main output)
# % required : no
# %end

# %option G_OPT_R_OUTPUT
# % key: circvarmap
# % description: Output Circular Variance layer
# % required : no
# %end

# %option G_OPT_R_OUTPUT
# % key: areamap
# % description: Output Area layer; values in square meters
# % required : no
# %end

# %option
# % key: thresh
# % type: double
# % description: Initial threshold (m^2).
# % required : yes
# %end

# %option
# % key: areamin
# % type: double
# % description: Minimum area (m^2) below whitch the slope unit is not further segmented
# % required : yes
# %end

# %option
# % key: areamax
# % type: double
# % description: Maximum area (m^2) above which the slope unit is segmented irrespective of aspect
# % required : no
# %end

# %option
# % key: cvmin
# % type: double
# % description: Minimum value of the circular variance (0.0-1.0) below whitch the slope unit is not further segmented
# % required : yes
# %end

# %option
# % key: rf
# % type: integer
# % description: Factor used to iterativelly reduce initial threshold: newthresh=thresh-thresh/reductionfactor
# % required : yes
# %end

# %option
# % key: maxiteration
# % type: integer
# % description: maximum number of iteration to do before the procedure is in any case stopped
# % required : yes
# %end

# %option
# % key: cleansize
# % type: double
# % description: Slope Units size to be removed
# % required : no
# %end

# %flag
# % key: m
# % description: Perform quick cleaning of small-sized areas and stripes
# % guisection: flags
# %end

# %flag
# % key: n
# % description: Perform detailed cleaning of small-sized areas (slow)
# % guisection: flags
# %end

import sys
import grass.script as grass
from grass.script import core as grasscore
import atexit
import os


def cleanup():
    if grass.find_file("slu_r")["file"]:
        grass.message(("Removing temporary files"))
        grass.run_command(
            "g.remove",
            type="raster",
            pattern=(
                "slu_diversity*,slopeslutmp,aspectslutmp,maxacc_*,cvar_*,count_*,stddevslope_*,slu_r*,todo*"
            ),
            flags="f",
            quiet=True,
        )
        grass.run_command(
            "g.remove",
            type="raster",
            pattern=(
                "slu_clump,slu_count,aspect,seno,coseno,sumseno,sumcoseno,count,cvar"
            ),
            flags="f",
            quiet=True,
        )
    if grass.find_file("MASK")["file"]:
        # grass.run_command('r.mask', flags='r')
        grass.run_command(
            "g.remove", type="raster", name=("MASK"), flags="f", quiet=True
        )


def main():
    global thc
    dem = options["demmap"]
    if options["plainsmap"]:
        plains = options["plainsmap"]
    slumap = options["slumap"]
    circvarmap = options["circvarmap"]
    areamap = options["areamap"]
    th = float(options["thresh"])
    amin = float(options["areamin"])
    cvarmin = float(options["cvmin"])
    red = int(options["rf"])
    maxiter = int(options["maxiteration"])
    if options["cleansize"]:
        cleansize = int(options["cleansize"])
        if options["slumapclean"]:
            slumapclean = options["slumapclean"]
        else:
            grass.fatal("When cleansize is provided, slumapclean is mandatory.")
        if flags["m"] and flags["n"]:
            grass.fatal(
                "When cleansize is provided, only one between m and n can be specified."
            )
    # estimating values of parameters in cells units
    region = grasscore.region()
    nsres = region["nsres"]
    ewres = region["ewres"]
    thc = int(th / (nsres * ewres))
    aminc = int(amin / (nsres * ewres))
    if options["areamax"]:
        # in square meters
        amax = float(options["areamax"])
        # in cells units
        amaxc = int(amax / (nsres * ewres))

    # setting the mask on the DTM
    grass.run_command("r.mask", raster=dem, overwrite=True, quiet=True)

    # generating the aspect layer
    grass.run_command(
        "r.slope.aspect",
        elevation=dem,
        aspect="aspectslutmp",
        overwrite=True,
        quiet=True,
    )

    # printing something
    grass.message(("Initial threshold (cells) is : %s") % thc)
    grass.message(("Initial minimum area (cells) is : %s") % aminc)

    # calculating sin and cosin of the aspect layer
    exp = "$out = cos($a)"
    grass.mapcalc(exp, out="coseno", a="aspectslutmp", overwrite=True, quiet=True)
    exp = "$out = sin($a)"
    grass.mapcalc(exp, out="seno", a="aspectslutmp", overwrite=True, quiet=True)

    # setting counters and base layers for the next "while loop"
    counter = 0
    last_counter = counter - 1
    control = 1
    i = grass.raster_info(dem)
    control_lastrun = int(i["cells"])
    grass.mapcalc("$out = null()", out="slu_r_0", overwrite=True, quiet=True)
    grass.mapcalc("$out = null()", out="cvar_0", overwrite=True, quiet=True)
    grass.mapcalc("$out = null()", out="count_0", overwrite=True, quiet=True)
    grass.mapcalc("$out = 1", out="slu_r_todo", overwrite=True, quiet=True)
    grass.mapcalc("$out = 1", out="slu_r_todo_0", overwrite=True, quiet=True)

    # starting the loop. The loop stops when: there are no halfbasins extracted (control=0),
    # OR the number of allowed iteration is exceeded (counter>maxiter) OR the threshold
    # (cells) is greather/equal than the reduction factor (thc >= red) otherwise int(thc-thc/red)
    # remains equal to thc
    while control > 0 and counter < maxiter and thc >= red:

        # generating the half-basins
        grass.run_command(
            "r.watershed",
            elevation=dem,
            hbasin="slu_r_tmp",
            thresh=thc,
            flags="abs",
            overwrite=True,
            quiet=True,
        )
        if options["plainsmap"]:
            exp = "$out = if(isnull($a),$b,null())"
            grass.mapcalc(
                exp, out="slu_r", a=plains, b="slu_r_tmp", overwrite=True, quiet=True
            )
        else:
            grass.run_command("g.copy", rast=("slu_r_tmp,slu_r"), quiet=True)

        grass.run_command("r.mask", raster="slu_r_todo", overwrite=True, quiet=True)
        grass.run_command(
            "r.stats.zonal",
            base="slu_r",
            cover="coseno",
            method="count",
            output="count",
            overwrite=True,
            quiet=True,
        )
        grass.run_command(
            "r.stats.zonal",
            base="slu_r",
            cover="coseno",
            method="sum",
            output="sumcoseno",
            overwrite=True,
            quiet=True,
        )
        grass.run_command(
            "r.stats.zonal",
            base="slu_r",
            cover="seno",
            method="sum",
            output="sumseno",
            overwrite=True,
            quiet=True,
        )

        # creating, for each half-basin, the layer where the circular variance is stored (cvar).
        # Circular variance is 1-R/n, R is the magnitude of the vectorial sum of all the unit
        # vectors of the aspect layer in each polygon and n is the number of unit vectors (and
        # cells) involved in the sum
        exp = "$out = 1-((sqrt(($a)^2 + ($b)^2))/$c)"
        grass.mapcalc(
            exp,
            out="cvar",
            a="sumseno",
            b="sumcoseno",
            c="count",
            overwrite=True,
            quiet=True,
        )
        grass.run_command("r.mask", flags="r", quiet=True)

        # selecting half-basins where area is larger than the minimum area and the average
        # unit vector is smaller than the unit vector threshold
        if options["areamax"]:
            exp = "$out = if($a>$f || ($a>$b && $c>$d),$g,null())"
            grass.mapcalc(
                exp,
                out="slu_r_todo",
                a="count",
                b=aminc,
                c="cvar",
                d=cvarmin,
                g="slu_r",
                f=amaxc,
                overwrite=True,
                quiet=True,
            )
            grass.run_command(
                "g.copy",
                rast=("count,count_prova_%s") % counter,
                quiet=True,
                overwrite=True,
            )
            # exp = "$out = if($a>$b,1,null())"
            # grass.mapcalc(exp, out = "slu_r_large"+str(counter), a = "count", b = amaxc , overwrite=True, quiet=True)
        else:
            exp = "$out = if($a>$b && $c>$d,$g,null())"
            grass.mapcalc(
                exp,
                out="slu_r_todo",
                a="count",
                b=aminc,
                c="cvar",
                d=cvarmin,
                g="slu_r",
                overwrite=True,
                quiet=True,
            )

        # checking that there actually are half-basins with area greater than areamin
        # and circular variance greater than cvarmin. otherwise the loop exits
        s = grass.read_command("r.univar", flags="g", map="slu_r_todo", quiet=True)
        kv = grass.parse_key_val(s)
        # ivan
        if kv["n"]:
            # if kv.has_key("n"):
            # increasing counter
            last_counter = counter
            counter = counter + 1
            # patching the current half-basins, cvar and counter that were not selected
            # in the previous steps with those that come from the previous step of the loop
            grass.run_command(
                "g.copy", rast=("slu_r_todo,slu_r_todo_%s") % counter, quiet=True
            )
            grass.run_command(
                "r.mask", raster="slu_r_todo", flags="i", overwrite=True, quiet=True
            )
            grass.run_command(
                "r.patch",
                input=("slu_r_" + str(last_counter), "slu_r"),
                output="slu_r_" + str(counter),
                overwrite=True,
                quiet=True,
            )
            grass.run_command(
                "g.copy",
                rast=("slu_r_" + str(counter), "slu_r_prova_" + str(counter)),
                quiet=True,
                overwrite=True,
            )
            grass.run_command(
                "r.patch",
                input=("cvar_" + str(last_counter), "cvar"),
                output="cvar_" + str(counter),
                overwrite=True,
                quiet=True,
            )
            grass.run_command(
                "r.patch",
                input=("count_" + str(last_counter), "count"),
                output="count_" + str(counter),
                overwrite=True,
                quiet=True,
            )
            grass.run_command("r.mask", flags="r", quiet=True)

            # rejecting partition if average area of new half-basins is less than amin;
            # not effective on large areas, if areamax is present
            if counter > 0:
                if options["areamax"]:
                    if counter == 1:
                        grass.mapcalc(
                            "$out = 1", out="count_prova_", overwrite=True, quiet=True
                        )
                    exp = "$out = if($a>$b,1,null())"
                    grass.mapcalc(
                        exp,
                        out="slu_r_large" + str(counter),
                        a="count_prova_" + str(last_counter - 1),
                        b=amaxc,
                        overwrite=True,
                        quiet=True,
                    )
                    exp = "$out = if(isnull($a),$b,null())"
                    grass.mapcalc(
                        exp,
                        out="MASK",
                        a="slu_r_large" + str(counter),
                        b="slu_r_" + str(counter),
                        overwrite=True,
                        quiet=True,
                    )
                    grass.run_command(
                        "g.copy",
                        rast=("MASK", "mask" + str(counter)),
                        quiet=True,
                        overwrite=True,
                    )
                else:
                    grass.run_command(
                        "r.mask", raster="slu_r_" + str(counter), quiet=True
                    )

                z = grass.read_command(
                    "r.univar",
                    flags="g",
                    map="slu_r_todo_" + str(last_counter),
                    quiet=True,
                )
                kvz = grass.parse_key_val(z)

                #                grass.message(("Univar: %s") % kvz )
                # ivan
                # if kvz.has_key("n"):
                if kvz["n"]:
                    en = int(kvz["n"])
                    grass.message(("Univar: %s") % en)
                    if en > 0:
                        grass.run_command(
                            "r.statistics",
                            base="slu_r_todo_" + str(last_counter),
                            cover="slu_r",
                            method="diversity",
                            output="slu_diversity_" + str(counter),
                            overwrite=True,
                            quiet=True,
                        )
                        #                        grass.run_command('r.univar', map='slu_r')
                        #                        grass.run_command('r.univar', map='slu_r_todo_'+str(last_counter))
                        #                        grass.run_command('r.univar', map='slu_diversity_'+str(counter))
                        grass.run_command(
                            "r.stats.zonal",
                            base="slu_r_todo_" + str(last_counter),
                            cover="coseno",
                            method="count",
                            output="todocount_" + str(counter),
                            overwrite=True,
                            quiet=True,
                        )
                        exp = "$out = $d/$a"
                        grass.mapcalc(
                            exp,
                            out="slu_r_test_" + str(counter),
                            a="@slu_diversity_" + str(counter),
                            d="todocount_" + str(counter),
                            overwrite=True,
                            quiet=True,
                        )
                        exp = "$out = if($d<$e,$c,null())"
                        grass.mapcalc(
                            exp,
                            out="slu_r_corr_" + str(counter),
                            b="slu_r_" + str(counter),
                            c="slu_r_todo_" + str(last_counter),
                            d="slu_r_test_" + str(counter),
                            e=aminc,
                            overwrite=True,
                            quiet=True,
                        )
                        grass.run_command("r.mask", flags="r", quiet=True)
                        grass.run_command(
                            "r.patch",
                            input=(
                                "slu_r_corr_" + str(counter),
                                "slu_r_" + str(counter),
                            ),
                            output="slu_r_" + str(counter),
                            overwrite=True,
                            quiet=True,
                        )
                    else:
                        grass.run_command("r.mask", flags="r", quiet=True)
                else:
                    grass.run_command("r.mask", flags="r", quiet=True)

            control = int(kv["n"])
            thc = int(thc - thc / red)
            thhect = thc * nsres * ewres / 10000
            grass.message(("Threshold (hectars) is: %s") % thhect)
            grass.message(
                ("No. of cells to be still classified as SLU is: %s. Loop done: %s")
                % (control, counter)
            )
        else:
            # exit the loop
            grass.message(("Nothing to do, ready to write the outputs "))
            control = 0

    # depending on how the while loop is exited the slu_r_$counter may have some small holes. Here we fill them.
    exp = "$out = if(isnull($a),$b,$a)"
    grass.mapcalc(
        exp,
        out="slu_r_" + str(counter),
        a="slu_r_" + str(counter),
        b="slu_r",
        overwrite=True,
        quiet=True,
    )
    exp = "$out = $a"
    grass.mapcalc(
        exp,
        out="cvar_" + str(counter),
        a="cvar_" + str(last_counter),
        overwrite=True,
        quiet=True,
    )
    exp = "$out = $a"
    grass.mapcalc(
        exp,
        out="count_" + str(counter),
        a="count_" + str(last_counter),
        overwrite=True,
        quiet=True,
    )

    # preparing the outputs
    exp = "$out = $a"
    grass.mapcalc(
        exp, out="slumap_1", a="slu_r_" + str(counter), overwrite=True, quiet=True
    )
    # add areas where DEM exists, and SUs do not exist
    if options["plainsmap"]:
        exp = "$out = if(isnull($a),if(isnull($b),if(isnull($c),null(),1),null()),$a)"
        grass.mapcalc(
            exp,
            a="slumap_1",
            b=plains,
            c=dem,
            out="slumap_2",
            overwrite=True,
            quiet=True,
        )
    else:
        exp = "$out = if(isnull($a),if(isnull($c),null(),1),$a)"
        grass.mapcalc(
            exp, a="slumap_1", c=dem, out="slumap_2", overwrite=True, quiet=True
        )
    grass.run_command(
        "r.clump", input="slumap_2", output=slumap, overwrite=True, quiet=True
    )
    grass.run_command("g.remove", type="raster", name="slumap_1", flags="f", quiet=True)
    grass.run_command("g.remove", type="raster", name="slumap_2", flags="f", quiet=True)
    grass.run_command(
        "r.colors", map="slu_r_" + str(counter), color="random", quiet=True
    )

    if circvarmap:
        print(circvarmap)
        exp = "$out = $a"
        grass.mapcalc(
            exp, out=circvarmap, a="cvar_" + str(counter), overwrite=True, quiet=True
        )
    if areamap:
        print(areamap)
        exp = "$out = $a*$b*$c"
        grass.mapcalc(
            exp,
            out=areamap,
            a="count_" + str(counter),
            b=nsres,
            c=ewres,
            overwrite=True,
            quiet=True,
        )
    if options["cleansize"]:
        if not flags["n"]:
            if not flags["m"]:
                print(" -- we want QUICK cleaning of small-sized areas: METHOD 1 --")
            grass.run_command("r.mask", raster=dem, overwrite=True, quiet=True)
            areamap = "areamap"
            #
            grass.run_command(
                "r.clump", input=slumap, output="slu_clump", overwrite=True, quiet=True
            )
            grass.run_command(
                "r.stats.zonal",
                base="slu_clump",
                cover="slu_clump",
                method="count",
                output="slu_count",
                overwrite=True,
                quiet=True,
            )
            #
            exp = "$out = $a*$b*$c"
            grass.mapcalc(
                exp,
                out=areamap,
                a="slu_count",
                b=nsres,
                c=ewres,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = if($a>$b,$c,null())"
            grass.mapcalc(
                exp,
                out="slu_r_clean",
                a=areamap,
                b=cleansize,
                c=slumap,
                overwrite=True,
                quiet=True,
            )
            cleansize = cleansize / (nsres * ewres)
            growdist = int((10 * cleansize / 3.14) ** 0.5)
            grass.run_command(
                "r.grow",
                input="slu_r_clean",
                output="slu_r_grow",
                radius=growdist,
                overwrite=True,
                quiet=True,
            )

        if flags["m"]:
            print(" -- we want QUICK cleaning of small-sized areas: METHOD 2 --")
            input = "slu_r_grow"
            output = "slu_no_stripes"
            grass.run_command(
                "r.neighbors",
                input=input,
                output="slu_diversity",
                method="diversity",
                size=5,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = if($a==1,1,null())"
            grass.mapcalc(
                exp,
                out="slu_diversity_nobordi",
                a="slu_diversity",
                overwrite=True,
                quiet=True,
            )
            grass.run_command(
                "r.grow",
                input="slu_diversity_nobordi",
                output="slu_diversity_nobordi_grow",
                radius=1.01,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = if(isnull($a),null(),$b)"
            grass.mapcalc(
                exp,
                out="slu_finale_nobordi",
                a="slu_diversity_nobordi_grow",
                b=input,
                overwrite=True,
                quiet=True,
            )
            grass.run_command(
                "r.grow",
                input="slu_finale_nobordi",
                output=output,
                radius=1000,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = int($a)"
            grass.mapcalc(exp, out=input, a=output, overwrite=True, quiet=True)
            grass.run_command(
                "g.remove",
                type="raster",
                name="slu_diversity,slu_diversity_nobordi,slu_diversity_nobordi_grow,slu_finale_nobordi,slu_no_stripes",
                flags="f",
                quiet=True,
            )
        if flags["n"]:
            print(" -- we want DETAILED cleaning of small-sized areas: METHOD 3 --")
            grass.run_command(
                "r.to.vect",
                input=slumap,
                output="slu_v_grow",
                type="area",
                overwrite=True,
                quiet=True,
            )
            grass.run_command(
                "r.slope.aspect",
                elevation=dem,
                aspect="aspect_slu",
                overwrite=True,
                quiet=True,
            )
            os.system(
                "/home/alvioli/dem_globo_srtm/range/script04/slu_code/clean_method_3.sh slu_v_grow vect2 %s"
                % cleansize
            )
            # applying method 2 at the end
            grass.run_command(
                "v.to.rast",
                input="vect2",
                output="rast2",
                use="cat",
                overwrite=True,
                quiet=True,
            )
            input = "rast2"
            output = "slu_r_grow"
            grass.run_command(
                "r.neighbors",
                input=input,
                output="slu_diversity",
                method="diversity",
                size=5,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = if($a==1,1,null())"
            grass.mapcalc(
                exp,
                out="slu_diversity_nobordi",
                a="slu_diversity",
                overwrite=True,
                quiet=True,
            )
            grass.run_command(
                "r.grow",
                input="slu_diversity_nobordi",
                output="slu_diversity_nobordi_grow",
                radius=1.01,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = if(isnull($a),null(),$b)"
            grass.mapcalc(
                exp,
                out="slu_finale_nobordi",
                a="slu_diversity_nobordi_grow",
                b=input,
                overwrite=True,
                quiet=True,
            )
            grass.run_command(
                "r.grow",
                input="slu_finale_nobordi",
                output=output,
                radius=1000,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = int($a)"
            grass.mapcalc(exp, out=input, a=output, overwrite=True, quiet=True)
            grass.run_command(
                "g.remove",
                type="raster",
                name="slu_diversity,slu_diversity_nobordi,slu_diversity_nobordi_grow,slu_finale_nobordi,rast2",
                flags="f",
                quiet=True,
            )
            grass.run_command(
                "g.remove",
                type="vector",
                name="slu_v_grow,vect2",
                flags="f",
                quiet=True,
            )

        if options["plainsmap"]:
            exp = "$out = if(isnull($b),if(isnull($c),null(),int($a)),null())"
            grass.mapcalc(
                exp,
                out=slumapclean,
                a="slu_r_grow",
                b=plains,
                c=dem,
                overwrite=True,
                quiet=True,
            )
        else:
            exp = "$out = if(isnull($c),null(),int($a))"
            grass.mapcalc(
                exp, out=slumapclean, a="slu_r_grow", c=dem, overwrite=True, quiet=True
            )
        grass.run_command("r.colors", map=slumapclean, color="random", quiet=True)


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
