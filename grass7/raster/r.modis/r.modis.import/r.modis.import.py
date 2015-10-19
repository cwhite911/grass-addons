#!/usr/bin/env python
# -*- coding: utf-8 -*-

############################################################################
#
# MODULE:        r.modis.import
# AUTHOR(S):     Luca Delucchi
# PURPOSE:       r.modis.import is an interface to pyModis for import into
#                GRASS GIS level 3 MODIS produts
#
# COPYRIGHT:        (C) 2011 by Luca Delucchi
#
#                This program is free software under the GNU General Public
#                License (>=v2). Read the file COPYING that comes with GRASS
#                for details.
#
#############################################################################
#
# REQUIREMENTS:
#   -  MRT tools, https://lpdaac.usgs.gov/lpdaac/tools/modis_reprojection_tool
#
#############################################################################

#%module
#% description: Import single or multiple tiles of MODIS products using pyModis/MRT
#% keyword: raster
#% keyword: MODIS
#%end
#%flag
#% key: m
#% description: Create a mosaic for each day
#%end
#%flag
#% key: t
#% description: Preserve temporary files (TIF and HDF mosaic)
#%end
#%flag
#% key: q
#% description: Ignore the QA map layer
#%end
#%option
#% key: mrtpath
#% type: string
#% key_desc: path
#% description: Full path to MRT directory
#% gisprompt: old,dir,input
#% required: no
#%end
#%option
#% key: dns
#% type: string
#% key_desc: path
#% description: Full path to single HDF file
#% gisprompt: old,file,input
#% required: no
#%end
#%option
#% key: files
#% type: string
#% key_desc: file
#% description: Full path to file with list of HDF files
#% gisprompt: old,file,input
#% required: no
#%end
#%option G_OPT_R_OUTPUT
#% required : no
#% guisection: Import
#%end
#%option
#% key: method
#% type: string
#% key_desc: resampling
#% description: Code of spatial resampling method
#% options: nearest, bilinear, cubic
#% answer: nearest
#% required: no
#%end
#%option
#% key: spectral
#% type: string
#% key_desc: spectral subset
#% description: String to choose subset of layers in HDF file for import
#% required: no
#%end

import os
import sys
import string
import glob
import shutil
import grass.script as grass
from datetime import date
from grass.pygrass.utils import get_lib_path

path = get_lib_path(modname='r.modis', libname='libmodis')
if path is None:
    grass.fatal("Not able to find the modis library directory.")
sys.path.append(path)

# try to import pymodis (modis) and some classes for r.modis.download
from rmodislib import resampling, product, projection, get_proj
from convertmodis import convertModis, createMosaic
from convertmodis_gdal import createMosaicGDAL, convertModisGDAL
from parsemodis import parseModis


def list_files(opt, mosaik=False):
    """If used in function single(): Return a list of HDF files from the file
    list. If used in function mosaic(): Return a dictionary with a list of HDF
    files for each day
    """
    # read the file with the list of HDF
    if opt['files'] != '':
        listoffile = open(opt['files'], 'r')
        basedir = os.path.split(listoffile.name)[0]
        # if mosaic create a dictionary
        if mosaik:
            filelist = {}
        # if not mosaic create a list
        else:
            filelist = []
        # append hdf files
        for line in listoffile:
            if string.find(line, 'xml') == -1 and mosaik is False:
                filelist.append(line.strip())
            # for mosaic create a list of hdf files for each day
            elif string.find(line, 'xml') == -1 and mosaik is True:
                day = line.split('/')[-1].split('.')[1]
                if day in filelist:
                    filelist[day].append(line)
                else:
                    filelist[day] = [line]
    # create a list for each file
    elif options['dns'] != '':
        filelist = [options['dns']]
        basedir = os.path.split(filelist[0])[0]
    return filelist, basedir


def spectral(opts, prod, q, m=False):
    """Return spectral string"""
    # return the spectral set selected by the user
    if opts['spectral'] != '':
        spectr = opts['spectral']
    # return the spectral by default
    else:
        if q:
            if prod['spec_qa']:
                spectr = prod['spec_qa']
            else:
                spectr = prod['spec']
        else:
            spectr = prod['spec']
        if m:
            spectr = spectr.replace(' 0', '')
    return spectr


def confile(pm, opts, q, mosaik=False):
    """Create the configuration file for MRT software"""
    # return projection and datum
    projObj = projection()
    proj = projObj.returned()
    dat = projObj.datum()
    if proj == 'UTM':
        zone = projObj.utmzone()
    else:
        zone = None
    cod = os.path.split(pm.hdfname)[1].split('.')[0]
    prod = product().fromcode(cod)
    if mosaik:
        # if mosaic it remove all the 0 from the subset string to convert all
        # the right layer
        spectr = spectral(opts, prod, q, True)
    else:
        spectr = spectral(opts, prod, q)
    # out prefix
    pref = prefix(opts)
    # resampling
    resampl = resampling(opts['method']).returned()
    # projpar
    projpar = projObj.return_params()
    if projpar != "( 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )":
        dat = "NODATUM"
    # resolution
    if proj != 'GEO':
        res = int(prod['res']) * int(projObj.proj['meters'])
    else:
        res = None
    try:
        conf = pm.confResample(spectr, res, pref, dat, resampl,
                               proj, zone, projpar)
        return conf
    except IOError, e:
        grass.fatal(e)


def prefix(options, name=False):
    """Return the prefix of output file if not set return None to use default
       value
    """
    if options['output']:
        return '%s.tif' % options['output']
    else:
        return None


def metadata(pars, mapp):
    """ Set metadata to the imported files """
    # metadata
    grass.run_command('r.support', quiet=True, map=mapp, source1="MODIS NASA",
                      hist="Imported with r.modis.import")
    # timestamp
    rangetime = pars.retRangeTime()
    data = rangetime['RangeBeginningDate'].split('-')
    dataobj = date(int(data[0]), int(data[1]), int(data[2]))
    grass.run_command('r.timestamp', map=mapp, quiet=True,
                      date=dataobj.strftime("%d %b %Y"))
    return 0
    # color
#    if string.find(mapp, 'QC') != -1 or string.find(mapp, 'Quality') != -1 or \
#    string.find(mapp, 'QA') != -1:
#        grass.run_command('r.colors', quiet=True, map=mapp, color=coll)
#    elif string.find(mapp, 'NDVI') != -1:
#        grass.run_command('r.colors', quiet=True, map=mapp, color=coll[0])
#    elif string.find(mapp, 'EVI') != -1:
#        grass.run_command('r.colors', quiet=True, map=mapp, color=coll[1])
#    elif string.find(mapp, 'LST') != -1:
#        grass.run_command('r.colors', quiet=True, map=mapp, color=coll[0])
#    elif string.find(mapp, 'Snow') != -1:
#        grass.run_command('r.colors', quiet=True, map=mapp, color=coll[0])


def import_tif(out, basedir, rem, write, pm, target=None):
    """Import TIF files"""
    print "start import"
    # list of tif files
    tifiles = glob.glob1(basedir, "*.tif")
    if not tifiles:
        tifiles = glob.glob1(os.getcwd(), "*.tif")
    if not tifiles:
        grass.fatal(_('Error during the conversion'))
    # check if user is in latlong location to set flag l
    if projection().val == 'll':
        f = "l"
    else:
        f = None
    outfile = []
    # for each file import it
    for t in tifiles:
        basename = os.path.splitext(t)[0]
        name = os.path.join(basedir, t)
        if not os.path.exists(name):
            name = os.path.join(os.getcwd(), t)
        if not os.path.exists(name):
            grass.warning(_("File %s doesn't find" % basename))
            continue
        filesize = int(os.path.getsize(name))
        if filesize < 1000:
            grass.warning(_('Probably some error occur during the conversion'
                            'for file <%s>. Escape import' % name))
            continue
        try:
            grass.run_command('r.in.gdal', input=name, output=basename,
                              flags=f, overwrite=write, quiet=True)
            outfile.append(basename)
        except:
            grass.warning(_('Error during import of %s' % basename))
            continue
        metadata(pm, basename)
        if rem:
            os.remove(name)
        if target:
            if target != basedir:
                shutil.move(name, target)
    return outfile


def findfile(pref, suff):
    """ Check if a file exists on mapset """
    if grass.find_file(pref + suff)['file']:
        return grass.find_file(pref + suff)
    else:
        grass.warning(_("Raster map <%s> not found") % (pref + suff))


def single(options, remove, an, ow):
    """Convert the HDF file to TIF and import it
    """
    listfile, basedir = list_files(options)
    # for each file
    for i in listfile:
        if os.path.exists(i):
            hdf = i
        else:
            # the full path to hdf file
            hdf = os.path.join(basedir, i)
            if not os.path.exists(hdf):
                grass.warning(_("%s not found" % i))
                continue
        if options['mrtpath']:
            # create conf file fro mrt tools
            pm = parseModis(hdf)
            confname = confile(pm, options, an)
            # create convertModis class and convert it in tif file
            execmodis = convertModis(hdf, confname, options['mrtpath'])
        else:
            projwkt = get_proj('w')
            projObj = projection()
            pref = listfile[0].split('/')[-1]
            prod = product().fromcode(pref.split('.')[0])
            spectr = spectral(options, prod, an)
            if projObj.returned() != 'GEO':
                res = int(prod['res']) * int(projObj.proj['meters'])
            else:
                res = None
            prod = product().fromcode(pref.split('.')[0])
            outname = "%s.%s.mosaic" % (pref.split('.')[0], pref.split('.')[1])
            execmodis = convertModisGDAL(hdf, outname, spectr, res,
                                         wkt=projwkt)
        execmodis.run()
        output = prefix(options)
        if not output:
            output = os.path.split(hdf)[1].rstrip('.hdf')
        # import tif files
        import_tif(output, basedir, remove, ow, pm)
        os.remove(confname)


def mosaic(options, remove, an, ow):
    """Create a daily mosaic of HDF files convert to TIF and import it
    """
    dictfile, targetdir = list_files(options, True)
    pid = str(os.getpid())
    # for each day
    for dat, listfiles in dictfile.iteritems():
        pref = listfiles[0].split('/')[-1]
        prod = product().fromcode(pref.split('.')[0])
        spectr = spectral(options, prod, an)
        spectr = spectr.lstrip('( ').rstrip(' )')
        outname = "%s.%s.mosaic" % (pref.split('.')[0], pref.split('.')[1])
        # create mosaic
        if options['mrtpath']:
            # create the file with the list of name
            tempfile = open(os.path.join(targetdir, pid), 'w')
            tempfile.writelines(listfiles)
            tempfile.close()
            # basedir of tempfile, where hdf files are write
            basedir = os.path.split(tempfile.name)[0]
            # return the spectral subset in according mrtmosaic tool format
            cm = createMosaic(tempfile.name, outname, options['mrtpath'],
                              spectr)
            cm.run()
        else:
            basedir = targetdir
            cm = createMosaicGDAL(listfiles, spectr)
            cm.write_vrt(outname)
        # list of hdf files
        hdfiles = glob.glob1(basedir, outname + "*.hdf")
        for i in hdfiles:
            # the full path to hdf file
            hdf = os.path.join(basedir, i)
            # create convertModis class and convert it in tif file
            if options['mrtpath']:
                # create conf file fro mrt tools
                pm = parseModis(hdf)
                confname = confile(pm, options, an, True)
                execmodis = convertModis(hdf, confname, options['mrtpath'])
            else:
                projwkt = get_proj('w')
                projObj = projection()
                if projwkt.returned() != 'GEO':
                    res = int(prod['res']) * int(projObj.proj['meters'])
                else:
                    res = None
                execmodis = convertModisGDAL(hdf, outname, res, wkt=projwkt,
                                             vrt=True)
            execmodis.run()
            # remove hdf
            if remove:
                # import tif files
                import_tif(outname, basedir, remove, ow, pm)
                os.remove(hdf)
                os.remove(hdf + '.xml')
            # move the hdf and hdf.xml to the dir where are the original files
            else:
                # import tif files
                import_tif(outname, basedir, remove, ow, pm, targetdir)
                try:
                    shutil.move(hdf, targetdir)
                    shutil.move(hdf + '.xml', targetdir)
                except:
                    pass
            # remove the conf file
            os.remove(confname)
        grass.try_remove(tempfile.name)
        grass.try_remove(os.path.join(targetdir, 'mosaic', pid))


def main():
    # check if you are in GRASS
    gisbase = os.getenv('GISBASE')
    if not gisbase:
        grass.fatal(_('$GISBASE not defined'))
        return 0
    # return an error if q and spectral are set
    if not flags['q'] and options['spectral'] != '':
        grass.warning(_('If no QA layer chosen in the "spectral" option'
                        ' the command will report an error'))
    # return an error if both dns and files option are set or not
    if options['dns'] == '' and options['files'] == '':
        grass.fatal(_('Choose one of "dns" or "files" options'))
        return 0
    elif options['dns'] != '' and options['files'] != '':
        grass.fatal(_('It is not possible set "dns" and "files"'
                      ' options together'))
        return 0
    # check the version
    version = grass.core.version()
    # this is would be set automatically
    if version['version'].find('7.') == -1:
        grass.fatal(_('GRASS GIS version 7 required'))
        return 0
    # check if remove the files or not
    if flags['t']:
        remove = False
    else:
        remove = True
    if grass.overwrite():
        over = True
    else:
        over = False
    # check if do check quality, rescaling and setting of colors
    if flags['q']:
        analyze = False
    else:
        analyze = True
    # check if import simple file or mosaic
    if flags['m'] and options['dns'] != '':
        grass.fatal(_('It is not possible to create a mosaic with a single'
                      ' HDF file'))
        return 0
    elif flags['m']:
        mosaic(options, remove, analyze, over)
    else:
        single(options, remove, analyze, over)

if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())
