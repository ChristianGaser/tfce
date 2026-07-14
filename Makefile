# Release the MATLAB/SPM toolbox.
#
# This Makefile is now only about the MATLAB side. The repository holds three
# things:
#
#   c/       the TFCE core, plain C, shared by both bindings and owned by neither
#   matlab/  the SPM toolbox: the m-files, the mex glue and the prebuilt binaries
#   python/  the pip package, released by .github/workflows/python_wheels.yml
#
# An installed SPM toolbox is a single flat folder, so the targets below assemble
# one: everything from matlab/ plus the headers from c/ that a user needs in order
# to recompile the mex-files. compile.m looks for the core beside itself first and
# in ../c second, which is exactly what makes both layouts work.
#
# $Id$

OLDVERSION="TFCE1.2"
NEWVERSION="TFCE1.3"
REVISION=`git rev-list --count HEAD`
DATE=`git log --date short |grep "Date:"|head -1|cut -f2 -d':'|sed -e s'/ //g'`
VERSION=`echo ${NEWVERSION} | sed -e 's/TFCE//g'`

TARGET=/Users/gaser/spm/spm12/toolbox/TFCE
TARGET2=/Volumes/UltraMax/spm12/toolbox/TFCE
TARGET3=paris.biomag.uni-jena.de:/home/gaser/spm12/toolbox/TFCE

STARGET_HOST=141.35.69.218
STARGET_FOLDER=/volume1/web/tfce
STARGET=${STARGET_HOST}:${STARGET_FOLDER}

# What goes into an installed toolbox. The mex sources and the core headers are
# shipped as well, so that a user on a platform we have no binary for can run
# compile themselves. tfce_batch.h is one of them: leave it out and the shipped
# sources no longer compile.
MATLAB_FILES=matlab/Contents.* matlab/tfce_*.m matlab/spm_TFCE.m matlab/snpm_P_FDR.m \
             matlab/tbx_cfg_tfce.m matlab/cat_spm_results_ui.m matlab/compile.m
MEX_FILES=matlab/tfceMex_maxtree.* matlab/tfceMex_maxtree_batch.* matlab/tfceMex_resss.*
CORE_FILES=c/tfce_maxtree.h c/tfce_batch.h c/tfce_threads.h
MISC_FILES=matlab/html

FILES=${MATLAB_FILES} ${MEX_FILES} ${CORE_FILES} ${MISC_FILES}

ZIPFILE=tfce_v${VERSION}.zip
ZIPFILE_OLD=tfce_r${REVISION}.zip

help:
	-@echo Available commands:
	-@echo "  install install2 install3   copy the toolbox into an SPM tree"
	-@echo "  doc                         regenerate matlab/html/tfce.html"
	-@echo "  update                      stamp version/revision into the m-files"
	-@echo "  zip                         build ${ZIPFILE}"
	-@echo "  scp                         tag, and upload the zip"
	-@echo ""
	-@echo "The Python package is released by .github/workflows/python_wheels.yml"
	-@echo "on a GitHub release, not from here."

install:
	-@echo install
	-@test ! -d ${TARGET} || rm -rf ${TARGET}
	-@mkdir -p ${TARGET}
	-@cp -R ${FILES} ${TARGET}

install2:
	-@echo install2
	-@test ! -d ${TARGET2} || rm -rf ${TARGET2}
	-@mkdir -p ${TARGET2}
	-@cp -R ${FILES} ${TARGET2}

install3:
	-@echo install3
	-@scp -r ${FILES} ${TARGET3}/

doc:
	-@cat matlab/html/tfce.txt | sed -e 's/VERSION/'${NEWVERSION}'/g' -e 's/RELNUMBER/r'${REVISION}'/g' -e 's/DATE/'${DATE}'/g' > matlab/html/tfce.html
	-@sed -i '' -e 's/VERSION/'${NEWVERSION}'/g' matlab/spm_TFCE.m

update: doc
	-@git fetch
	-@echo '% TFCE Toolbox' > matlab/Contents.m
	-@echo '% Version ' ${REVISION} ' ('${NEWVERSION}')' ${DATE} >> matlab/Contents.m
	-@cat matlab/Contents_info.txt >> matlab/Contents.m
	-@perl -p -i -e "s/${OLDVERSION}/${NEWVERSION}/g" matlab/spm_TFCE.m
	-@echo '% TFCE Toolbox' > matlab/INSTALL.txt
	-@echo '% Version ' ${REVISION} ${NEWVERSION} ${DATE} >> matlab/INSTALL.txt
	-@cat matlab/INSTALL_info.txt >> matlab/INSTALL.txt

# The zip is the flat folder SPM expects: matlab/ and c/ collapse into one TFCE/.
zip: update
	-@echo zip
	-@test ! -d TFCE || rm -r TFCE
	-@mkdir TFCE
	-@cp -rp ${FILES} TFCE
	-@bash update_revision.sh
	-@zip ${ZIPFILE} -rm TFCE

scp: zip
	-@git tag -f ${VERSION} -m "Release ${VERSION}"
	-@echo scp to http://${STARGET_HOST}/tfce/${ZIPFILE}
	-@scp -O -P 2222 CHANGES.txt ${ZIPFILE} ${STARGET}/${ZIPFILE_OLD}
	-@bash -c "ssh -p 2222 ${STARGET_HOST} ln -fs ${STARGET_FOLDER}/${ZIPFILE_OLD} ${STARGET_FOLDER}/tfce_latest.zip"

.PHONY: help install install2 install3 doc update zip scp
