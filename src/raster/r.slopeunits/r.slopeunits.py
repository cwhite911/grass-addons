#!/usr/bin/env python3
#
############################################################################
#
# MODULE:	r.slopeunits for GRASS 8
# AUTHOR(S):    Ivan Marchesini, Massimiliano Alvioli, Corey T. White (Implemented in GRASS Addons)
# PURPOSE:	To create a raster layer of slope units
# COPYRIGHT: (C) 2004-2024 by the GRASS Development Team
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
import grass.script as gs
from grass.script import core as grasscore
import grass.pygrass.vector as vector
import grass.pygrass.raster as raster
from grass.script import mapcalc
import atexit
import os


def cleanup():
    if gs.find_file("slu_r")["file"]:
        gs.message(("Removing temporary files"))
        gs.run_command(
            "g.remove",
            type="raster",
            pattern=(
                "slu_diversity*,slopeslutmp,aspectslutmp,maxacc_*,cvar_*,count_*,stddevslope_*,slu_r*,todo*"
            ),
            flags="f",
            quiet=True,
        )
        gs.run_command(
            "g.remove",
            type="raster",
            pattern=(
                "slu_clump,slu_count,aspect,seno,coseno,sumseno,sumcoseno,count,cvar"
            ),
            flags="f",
            quiet=True,
        )
    if gs.find_file("MASK")["file"]:
        # gs.run_command('r.mask', flags='r')
        gs.run_command("g.remove", type="raster", name=("MASK"), flags="f", quiet=True)


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
            gs.fatal("When cleansize is provided, slumapclean is mandatory.")
        if flags["m"] and flags["n"]:
            gs.fatal(
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
    gs.run_command("r.mask", raster=dem, overwrite=True, quiet=True)

    # generating the aspect layer
    gs.run_command(
        "r.slope.aspect",
        elevation=dem,
        aspect="aspectslutmp",
        overwrite=True,
        quiet=True,
    )

    # printing something
    gs.message(("Initial threshold (cells) is : %s") % thc)
    gs.message(("Initial minimum area (cells) is : %s") % aminc)

    # calculating sin and cosin of the aspect layer
    exp = "$out = cos($a)"
    gs.mapcalc(exp, out="coseno", a="aspectslutmp", overwrite=True, quiet=True)
    exp = "$out = sin($a)"
    gs.mapcalc(exp, out="seno", a="aspectslutmp", overwrite=True, quiet=True)

    # setting counters and base layers for the next "while loop"
    counter = 0
    last_counter = counter - 1
    control = 1
    i = gs.raster_info(dem)
    control_lastrun = int(i["cells"])
    gs.mapcalc("$out = null()", out="slu_r_0", overwrite=True, quiet=True)
    gs.mapcalc("$out = null()", out="cvar_0", overwrite=True, quiet=True)
    gs.mapcalc("$out = null()", out="count_0", overwrite=True, quiet=True)
    gs.mapcalc("$out = 1", out="slu_r_todo", overwrite=True, quiet=True)
    gs.mapcalc("$out = 1", out="slu_r_todo_0", overwrite=True, quiet=True)

    # starting the loop. The loop stops when: there are no halfbasins extracted (control=0),
    # OR the number of allowed iteration is exceeded (counter>maxiter) OR the threshold
    # (cells) is greather/equal than the reduction factor (thc >= red) otherwise int(thc-thc/red)
    # remains equal to thc
    while control > 0 and counter < maxiter and thc >= red:

        # generating the half-basins
        gs.run_command(
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
            gs.mapcalc(
                exp, out="slu_r", a=plains, b="slu_r_tmp", overwrite=True, quiet=True
            )
        else:
            gs.run_command("g.copy", rast=("slu_r_tmp,slu_r"), quiet=True)

        gs.run_command("r.mask", raster="slu_r_todo", overwrite=True, quiet=True)
        gs.run_command(
            "r.stats.zonal",
            base="slu_r",
            cover="coseno",
            method="count",
            output="count",
            overwrite=True,
            quiet=True,
        )
        gs.run_command(
            "r.stats.zonal",
            base="slu_r",
            cover="coseno",
            method="sum",
            output="sumcoseno",
            overwrite=True,
            quiet=True,
        )
        gs.run_command(
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
        gs.mapcalc(
            exp,
            out="cvar",
            a="sumseno",
            b="sumcoseno",
            c="count",
            overwrite=True,
            quiet=True,
        )
        gs.run_command("r.mask", flags="r", quiet=True)

        # selecting half-basins where area is larger than the minimum area and the average
        # unit vector is smaller than the unit vector threshold
        if options["areamax"]:
            exp = "$out = if($a>$f || ($a>$b && $c>$d),$g,null())"
            gs.mapcalc(
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
            gs.run_command(
                "g.copy",
                rast=("count,count_prova_%s") % counter,
                quiet=True,
                overwrite=True,
            )
            # exp = "$out = if($a>$b,1,null())"
            # gs.mapcalc(exp, out = "slu_r_large"+str(counter), a = "count", b = amaxc , overwrite=True, quiet=True)
        else:
            exp = "$out = if($a>$b && $c>$d,$g,null())"
            gs.mapcalc(
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
        s = gs.read_command("r.univar", flags="g", map="slu_r_todo", quiet=True)
        kv = gs.parse_key_val(s)
        # ivan
        if kv["n"]:
            # if kv.has_key("n"):
            # increasing counter
            last_counter = counter
            counter = counter + 1
            # patching the current half-basins, cvar and counter that were not selected
            # in the previous steps with those that come from the previous step of the loop
            gs.run_command(
                "g.copy", rast=("slu_r_todo,slu_r_todo_%s") % counter, quiet=True
            )
            gs.run_command(
                "r.mask", raster="slu_r_todo", flags="i", overwrite=True, quiet=True
            )
            gs.run_command(
                "r.patch",
                input=("slu_r_" + str(last_counter), "slu_r"),
                output="slu_r_" + str(counter),
                overwrite=True,
                quiet=True,
            )
            gs.run_command(
                "g.copy",
                rast=("slu_r_" + str(counter), "slu_r_prova_" + str(counter)),
                quiet=True,
                overwrite=True,
            )
            gs.run_command(
                "r.patch",
                input=("cvar_" + str(last_counter), "cvar"),
                output="cvar_" + str(counter),
                overwrite=True,
                quiet=True,
            )
            gs.run_command(
                "r.patch",
                input=("count_" + str(last_counter), "count"),
                output="count_" + str(counter),
                overwrite=True,
                quiet=True,
            )
            gs.run_command("r.mask", flags="r", quiet=True)

            # rejecting partition if average area of new half-basins is less than amin;
            # not effective on large areas, if areamax is present
            if counter > 0:
                if options["areamax"]:
                    if counter == 1:
                        gs.mapcalc(
                            "$out = 1", out="count_prova_", overwrite=True, quiet=True
                        )
                    exp = "$out = if($a>$b,1,null())"
                    gs.mapcalc(
                        exp,
                        out="slu_r_large" + str(counter),
                        a="count_prova_" + str(last_counter - 1),
                        b=amaxc,
                        overwrite=True,
                        quiet=True,
                    )
                    exp = "$out = if(isnull($a),$b,null())"
                    gs.mapcalc(
                        exp,
                        out="MASK",
                        a="slu_r_large" + str(counter),
                        b="slu_r_" + str(counter),
                        overwrite=True,
                        quiet=True,
                    )
                    gs.run_command(
                        "g.copy",
                        rast=("MASK", "mask" + str(counter)),
                        quiet=True,
                        overwrite=True,
                    )
                else:
                    gs.run_command("r.mask", raster="slu_r_" + str(counter), quiet=True)

                z = gs.read_command(
                    "r.univar",
                    flags="g",
                    map="slu_r_todo_" + str(last_counter),
                    quiet=True,
                )
                kvz = gs.parse_key_val(z)

                #                gs.message(("Univar: %s") % kvz )
                # ivan
                # if kvz.has_key("n"):
                if kvz["n"]:
                    en = int(kvz["n"])
                    gs.message(("Univar: %s") % en)
                    if en > 0:
                        gs.run_command(
                            "r.statistics",
                            base="slu_r_todo_" + str(last_counter),
                            cover="slu_r",
                            method="diversity",
                            output="slu_diversity_" + str(counter),
                            overwrite=True,
                            quiet=True,
                        )
                        #                        gs.run_command('r.univar', map='slu_r')
                        #                        gs.run_command('r.univar', map='slu_r_todo_'+str(last_counter))
                        #                        gs.run_command('r.univar', map='slu_diversity_'+str(counter))
                        gs.run_command(
                            "r.stats.zonal",
                            base="slu_r_todo_" + str(last_counter),
                            cover="coseno",
                            method="count",
                            output="todocount_" + str(counter),
                            overwrite=True,
                            quiet=True,
                        )
                        exp = "$out = $d/$a"
                        gs.mapcalc(
                            exp,
                            out="slu_r_test_" + str(counter),
                            a="@slu_diversity_" + str(counter),
                            d="todocount_" + str(counter),
                            overwrite=True,
                            quiet=True,
                        )
                        exp = "$out = if($d<$e,$c,null())"
                        gs.mapcalc(
                            exp,
                            out="slu_r_corr_" + str(counter),
                            b="slu_r_" + str(counter),
                            c="slu_r_todo_" + str(last_counter),
                            d="slu_r_test_" + str(counter),
                            e=aminc,
                            overwrite=True,
                            quiet=True,
                        )
                        gs.run_command("r.mask", flags="r", quiet=True)
                        gs.run_command(
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
                        gs.run_command("r.mask", flags="r", quiet=True)
                else:
                    gs.run_command("r.mask", flags="r", quiet=True)

            control = int(kv["n"])
            thc = int(thc - thc / red)
            thhect = thc * nsres * ewres / 10000
            gs.message(("Threshold (hectars) is: %s") % thhect)
            gs.message(
                ("No. of cells to be still classified as SLU is: %s. Loop done: %s")
                % (control, counter)
            )
        else:
            # exit the loop
            gs.message(("Nothing to do, ready to write the outputs "))
            control = 0

    # depending on how the while loop is exited the slu_r_$counter may have some small holes. Here we fill them.
    exp = "$out = if(isnull($a),$b,$a)"
    gs.mapcalc(
        exp,
        out="slu_r_" + str(counter),
        a="slu_r_" + str(counter),
        b="slu_r",
        overwrite=True,
        quiet=True,
    )
    exp = "$out = $a"
    gs.mapcalc(
        exp,
        out="cvar_" + str(counter),
        a="cvar_" + str(last_counter),
        overwrite=True,
        quiet=True,
    )
    exp = "$out = $a"
    gs.mapcalc(
        exp,
        out="count_" + str(counter),
        a="count_" + str(last_counter),
        overwrite=True,
        quiet=True,
    )

    # preparing the outputs
    exp = "$out = $a"
    gs.mapcalc(
        exp, out="slumap_1", a="slu_r_" + str(counter), overwrite=True, quiet=True
    )
    # add areas where DEM exists, and SUs do not exist
    if options["plainsmap"]:
        exp = "$out = if(isnull($a),if(isnull($b),if(isnull($c),null(),1),null()),$a)"
        gs.mapcalc(
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
        gs.mapcalc(exp, a="slumap_1", c=dem, out="slumap_2", overwrite=True, quiet=True)
    gs.run_command(
        "r.clump", input="slumap_2", output=slumap, overwrite=True, quiet=True
    )
    gs.run_command("g.remove", type="raster", name="slumap_1", flags="f", quiet=True)
    gs.run_command("g.remove", type="raster", name="slumap_2", flags="f", quiet=True)
    gs.run_command("r.colors", map="slu_r_" + str(counter), color="random", quiet=True)

    if circvarmap:
        print(circvarmap)
        exp = "$out = $a"
        gs.mapcalc(
            exp, out=circvarmap, a="cvar_" + str(counter), overwrite=True, quiet=True
        )
    if areamap:
        print(areamap)
        exp = "$out = $a*$b*$c"
        gs.mapcalc(
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
            gs.run_command("r.mask", raster=dem, overwrite=True, quiet=True)
            areamap = "areamap"
            #
            gs.run_command(
                "r.clump", input=slumap, output="slu_clump", overwrite=True, quiet=True
            )
            gs.run_command(
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
            gs.mapcalc(
                exp,
                out=areamap,
                a="slu_count",
                b=nsres,
                c=ewres,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = if($a>$b,$c,null())"
            gs.mapcalc(
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
            gs.run_command(
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
            gs.run_command(
                "r.neighbors",
                input=input,
                output="slu_diversity",
                method="diversity",
                size=5,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = if($a==1,1,null())"
            gs.mapcalc(
                exp,
                out="slu_diversity_nobordi",
                a="slu_diversity",
                overwrite=True,
                quiet=True,
            )
            gs.run_command(
                "r.grow",
                input="slu_diversity_nobordi",
                output="slu_diversity_nobordi_grow",
                radius=1.01,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = if(isnull($a),null(),$b)"
            gs.mapcalc(
                exp,
                out="slu_finale_nobordi",
                a="slu_diversity_nobordi_grow",
                b=input,
                overwrite=True,
                quiet=True,
            )
            gs.run_command(
                "r.grow",
                input="slu_finale_nobordi",
                output=output,
                radius=1000,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = int($a)"
            gs.mapcalc(exp, out=input, a=output, overwrite=True, quiet=True)
            gs.run_command(
                "g.remove",
                type="raster",
                name="slu_diversity,slu_diversity_nobordi,slu_diversity_nobordi_grow,slu_finale_nobordi",
                flags="f",
                quiet=True,
            )
        if flags["n"]:
            print(" -- we want DETAILED cleaning of small-sized areas: METHOD 3 --")
            gs.run_command(
                "r.to.vect",
                input=slumap,
                output="slu_v_grow",
                type="area",
                overwrite=True,
                quiet=True,
            )
            gs.run_command(
                "r.slope.aspect",
                elevation=dem,
                aspect="aspect_slu",
                overwrite=True,
                quiet=True,
            )
            # This was a shell called clean_method_3.sh
            clean_method_3("slu_v_grow", "vect2", cleansize)
            # applying method 2 at the end
            gs.run_command(
                "v.to.rast",
                input="vect2",
                output="rast2",
                use="cat",
                overwrite=True,
                quiet=True,
            )
            input = "rast2"
            output = "slu_r_grow"
            gs.run_command(
                "r.neighbors",
                input=input,
                output="slu_diversity",
                method="diversity",
                size=5,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = if($a==1,1,null())"
            gs.mapcalc(
                exp,
                out="slu_diversity_nobordi",
                a="slu_diversity",
                overwrite=True,
                quiet=True,
            )
            gs.run_command(
                "r.grow",
                input="slu_diversity_nobordi",
                output="slu_diversity_nobordi_grow",
                radius=1.01,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = if(isnull($a),null(),$b)"
            gs.mapcalc(
                exp,
                out="slu_finale_nobordi",
                a="slu_diversity_nobordi_grow",
                b=input,
                overwrite=True,
                quiet=True,
            )
            gs.run_command(
                "r.grow",
                input="slu_finale_nobordi",
                output=output,
                radius=1000,
                overwrite=True,
                quiet=True,
            )
            exp = "$out = int($a)"
            gs.mapcalc(exp, out=input, a=output, overwrite=True, quiet=True)
            gs.run_command(
                "g.remove",
                type="raster",
                name="slu_diversity,slu_diversity_nobordi,slu_diversity_nobordi_grow,slu_finale_nobordi,rast2",
                flags="f",
                quiet=True,
            )
            gs.run_command(
                "g.remove",
                type="vector",
                name="slu_v_grow,vect2",
                flags="f",
                quiet=True,
            )

        if options["plainsmap"]:
            exp = "$out = if(isnull($b),if(isnull($c),null(),int($a)),null())"
            gs.mapcalc(
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
            gs.mapcalc(
                exp, out=slumapclean, a="slu_r_grow", c=dem, overwrite=True, quiet=True
            )
        gs.run_command("r.colors", map=slumapclean, color="random", quiet=True)


def clean_method_3(input_vect, output_vect, minarea):
    # Get the current region
    region = gs.parse_command("g.region", flags="pg")
    nsres = float(region["nsres"])
    ewres = float(region["ewres"])

    smarea = 10 * nsres * ewres
    gs.run_command(
        "v.clean",
        input=input_vect,
        output="slu_clean",
        tool="rmarea",
        threshold=smarea,
        overwrite=True,
        quiet=True,
    )

    gs.run_command(
        "v.db.addcolumn",
        map="slu_clean",
        columns="area double, perimetro double",
        quiet=True,
    )
    gs.run_command(
        "v.to.db", map="slu_clean", option="area", columns="area", quiet=True
    )
    gs.run_command(
        "v.to.db", map="slu_clean", option="perimeter", columns="perimetro", quiet=True
    )
    gs.run_command(
        "v.db.droprow",
        input="slu_clean",
        where='area=""',
        output="slu_area",
        overwrite=True,
        quiet=True,
    )

    # Select features below the minimum area
    lista = gs.read_command(
        "v.db.select",
        flags="c",
        map="slu_area",
        where=f"area<={minarea}",
        columns="cat",
    ).strip()
    buchi = ",".join(lista.split())
    totalebuchi = len(buchi.split(","))

    gs.run_command(
        "v.extract",
        input="slu_area",
        cats=buchi,
        output="slu_buchi",
        type="area",
        overwrite=True,
        quiet=True,
    )
    gs.run_command(
        "v.to.rast",
        input="slu_buchi",
        output="slu_buchi",
        use="cat",
        overwrite=True,
        quiet=True,
    )

    # Select features above the minimum area
    lista = gs.read_command(
        "v.db.select", flags="c", map="slu_area", where=f"area>{minarea}", columns="cat"
    ).strip()
    nobuchi = ",".join(lista.split())

    # TODO: What is slu_nobuchi?
    # gs.run_command('v.to.rast', input='slu_nobuchi', output='slu_nobuchi', use='cat', overwrite=True, quiet=True)
    # gs.run_command('v.extract', input='slu_area', cats=nobuchi, output='slu_nobuchi', type='area', overwrite=True, quiet=True)

    gs.run_command(
        "v.category",
        input="slu_area",
        output="slu_bordi",
        layer=2,
        type="boundary",
        option="add",
        overwrite=True,
        quiet=True,
    )
    gs.run_command(
        "v.db.addtable",
        map="slu_bordi",
        layer=2,
        columns="left integer, right integer, lunghezza double",
        quiet=True,
    )
    gs.run_command(
        "v.to.db",
        map="slu_bordi",
        option="sides",
        columns="left,right",
        layer=2,
        type="boundary",
        quiet=True,
    )
    gs.run_command(
        "v.to.db",
        map="slu_bordi",
        option="length",
        columns="lunghezza",
        layer=2,
        type="boundary",
        quiet=True,
    )
    gs.run_command(
        "v.to.rast",
        input="slu_area",
        output="slu_area",
        use="cat",
        overwrite=True,
        quiet=True,
    )

    mapcalc("coseno = cos(aspect_slu)", overwrite=True, quiet=True)
    mapcalc("seno = sin(aspect_slu)", overwrite=True, quiet=True)

    # r.stats.zonal replacements for gs 7
    gs.run_command(
        "r.stats.zonal",
        base="slu_area",
        cover="coseno",
        method="count",
        output="count",
        overwrite=True,
        quiet=True,
    )
    gs.run_command(
        "r.stats.zonal",
        base="slu_area",
        cover="coseno",
        method="sum",
        output="sumcos",
        overwrite=True,
        quiet=True,
    )
    gs.run_command(
        "r.stats.zonal",
        base="slu_area",
        cover="seno",
        method="sum",
        output="sumsin",
        overwrite=True,
        quiet=True,
    )

    mapcalc("cos_medio = sumcos / count", overwrite=True, quiet=True)
    mapcalc("sin_medio = sumsin / count", overwrite=True, quiet=True)

    gs.run_command(
        "r.to.vect",
        input="cos_medio",
        output="cos_medio",
        type="area",
        overwrite=True,
        quiet=True,
    )
    gs.run_command(
        "r.to.vect",
        input="sin_medio",
        output="sin_medio",
        type="area",
        overwrite=True,
        quiet=True,
    )

    gs.run_command(
        "v.overlay",
        ainput="slu_area",
        binput="cos_medio",
        operator="and",
        atype="area",
        btype="area",
        output="cos_medio_over",
        overwrite=True,
        quiet=True,
    )
    gs.run_command(
        "v.overlay",
        ainput="slu_area",
        binput="sin_medio",
        operator="and",
        atype="area",
        btype="area",
        output="sin_medio_over",
        overwrite=True,
        quiet=True,
    )

    # Clean up areas
    pulire = (
        gs.read_command("v.category", input="slu_buchi", option="print", quiet=True)
        .strip()
        .split()
    )

    gs.run_command(
        "g.copy", vector=f"slu_area,{output_vect}", overwrite=True, quiet=True
    )

    ico = 1
    for i in pulire:
        lista1 = gs.read_command(
            "db.select",
            sql=f"select b2.right from slu_bordi_2 b2 where b2.left={i} and b2.right<>-1",
        ).strip()
        lista2 = gs.read_command(
            "db.select",
            sql=f"select b2.left from slu_bordi_2 b2 where b2.left<>-1 and b2.right={i}",
        ).strip()
        vicini = ",".join(sorted(set(lista1.split() + lista2.split()))[:-2])

        if vicini:
            print(
                f" --- --- -- buco numero {ico} di {totalebuchi}, cat: {i}, vicini: {vicini}"
            )
            ico += 1
            gs.run_command(
                "v.extract",
                input=output_vect,
                cats=vicini,
                output="intorno",
                type="area",
                overwrite=True,
                quiet=True,
            )
            chk_intorno = gs.read_command(
                "v.category",
                input="intorno",
                type="centroid",
                option="print",
                quiet=True,
            ).strip()

            if chk_intorno:
                gs.run_command(
                    "v.overlay",
                    ainput="intorno",
                    binput="slu_nobuchi",
                    output="intorno_OK",
                    atype="area",
                    btype="area",
                    olayer=[0, 1, 0],
                    operator="and",
                    overwrite=True,
                    quiet=True,
                )
                cos_buco = gs.read_command(
                    "v.db.select",
                    map="cos_medio_over",
                    where=f"a_cat={i}",
                    columns="b_value",
                    quiet=True,
                ).strip()
                sin_buco = gs.read_command(
                    "v.db.select",
                    map="sin_medio_over",
                    where=f"a_cat={i}",
                    columns="b_value",
                    quiet=True,
                ).strip()

                massimo = -10000
                jmax = 0
                loop = (
                    gs.read_command(
                        "v.category", input="intorno_OK", option="print", quiet=True
                    )
                    .strip()
                    .split()
                )

                for j in loop:
                    cos_j = gs.read_command(
                        "v.db.select",
                        map="cos_medio_over",
                        where=f"a_cat={j}",
                        columns="b_value",
                        quiet=True,
                    ).strip()
                    sin_j = gs.read_command(
                        "v.db.select",
                        map="sin_medio_over",
                        where=f"a_cat={j}",
                        columns="b_value",
                        quiet=True,
                    ).strip()
                    dotpr = (
                        int(
                            float(cos_buco) * float(cos_j)
                            + float(sin_buco) * float(sin_j)
                        )
                        * 10000
                    )

                    if dotpr >= massimo and dotpr > 0:
                        massimo = dotpr
                        jmax = j

                    print(i, j, cos_j, sin_j, dotpr, jmax)

                print(f"massimo: {massimo} per j={jmax}")

                if int(jmax) > 0:
                    lunghezza = gs.read_command(
                        "db.select",
                        sql=f"select b2.lunghezza from slu_bordi_2 b2 where (b2.left={i} and b2.right={jmax}) or (b2.left={jmax} and b2.right={i})",
                    ).strip()
                    perimetro = gs.read_command(
                        "v.db.select",
                        map="slu_clean",
                        columns="perimetro",
                        where=f"cat={i}",
                        quiet=True,
                    ).strip()

                    if int(lunghezza) > 0 and int(perimetro) > 0:
                        frazione = int(int(lunghezza) / int(perimetro) * 10000)

                        if frazione > 500:
                            print(
                                f"lungh: {lunghezza}; perim: {perimetro}; fract: {frazione}"
                            )
                            gs.run_command(
                                "v.extract",
                                input=output_vect,
                                output="slu_i",
                                cats=f"{i},{jmax}",
                                flags="d",
                                overwrite=True,
                                quiet=True,
                            )
                            gs.run_command(
                                "v.overlay",
                                ainput=output_vect,
                                binput="slu_i",
                                atype="area",
                                btype="area",
                                operator="not",
                                olayer=[0, 1, 0],
                                output="slu_j",
                                overwrite=True,
                                quiet=True,
                            )
                            gs.run_command(
                                "v.overlay",
                                ainput="slu_i",
                                binput="slu_j",
                                atype="area",
                                btype="area",
                                operator="or",
                                output="slu_k",
                                olayer=[1, 0, 0],
                                overwrite=True,
                                quiet=True,
                            )

                            gs.run_command(
                                "v.db.addcolumn",
                                map="slu_k",
                                column="newcat integer",
                                overwrite=True,
                                quiet=True,
                            )
                            gs.run_command(
                                "v.db.update",
                                map="slu_k",
                                layer=1,
                                column="newcat",
                                qcolumn="a_cat",
                                where="a_cat is not null",
                                overwrite=True,
                                quiet=True,
                            )
                            gs.run_command(
                                "v.db.update",
                                map="slu_k",
                                layer=1,
                                column="newcat",
                                qcolumn="b_cat",
                                where="b_cat is not null",
                                overwrite=True,
                                quiet=True,
                            )
                            gs.run_command(
                                "v.reclass",
                                input="slu_k",
                                output=output_vect,
                                column="newcat",
                                overwrite=True,
                                quiet=True,
                            )
                            gs.run_command(
                                "v.db.addtable",
                                map=output_vect,
                                overwrite=True,
                                quiet=True,
                            )
                            gs.run_command(
                                "v.db.addcolumn",
                                map=output_vect,
                                columns="area integer",
                                overwrite=True,
                                quiet=True,
                            )
                            gs.run_command(
                                "v.to.db",
                                map=output_vect,
                                option="area",
                                columns="area",
                                overwrite=True,
                                quiet=True,
                            )
                            gs.run_command(
                                "g.remove",
                                type="vector",
                                name=["slu_i", "slu_j", "slu_k"],
                                flags="f",
                                quiet=True,
                            )

            gs.run_command(
                "g.remove",
                type="vector",
                name=["intorno", "intorno_OK"],
                flags="f",
                quiet=True,
            )


if __name__ == "__main__":
    options, flags = gs.parser()
    atexit.register(cleanup)
    main()
