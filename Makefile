# Personal Makefile variables
#
# $Id$

OLDVERSION="TFCE1.0"
NEWVERSION="TFCE1.1"
REVISION=`git rev-list --count HEAD`
DATE=`git log --date short |grep "Date:"|head -1|cut -f2 -d':'|sed -e s'/ //g'`
VERSION=`echo ${NEWVERSION} | sed -e 's/TFCE//g'`

ZIPFOLDER=/Users/gaser/matlab/tfce8

TARGET=/Users/gaser/spm/spm12/toolbox/TFCE
TARGET2=/Volumes/UltraMax/spm12/toolbox/TFCE
TARGET3=paris.biomag.uni-jena.de:/home/gaser/spm12/toolbox/TFCE

STARGET_HOST=141.35.69.218
STARGET_HTDOCS=${STARGET_HOST}:/volume1/web/
STARGET_FOLDER=/volume1/web/tfce

STARGET=${STARGET_HOST}:${STARGET_FOLDER}

MATLAB_FILES=Contents.* tfce_*.m spm_TFCE.m snpm_P_FDR.m tbx_cfg_tfce.m cat_spm_results_ui.m
C_FILES=tfceMex_pthread.* 
MISC_FILES=html

FILES=${MATLAB_FILES} ${C_FILES} ${MISC_FILES}

ZIPFILE=tfce_r${REVISION}.zip

install:
	-@echo install
	-@test ! -d ${TARGET} || rm -rf ${TARGET}
	-@mkdir ${TARGET}
	-@cp -R ${FILES} ${TARGET}

install2:
	-@echo install2
	-@test ! -d ${TARGET2} || rm -rf ${TARGET2}
	-@mkdir ${TARGET2}
	-@cp -R ${FILES} ${TARGET2}

install3:
	-@echo install3
	-@scp -r ${FILES} ${TARGET3}/

help:
	-@echo Available commands:
	-@echo install install2 install3 zip scp doc update

doc:
	-@cat html/tfce.txt | sed -e 's/VERSION/'${NEWVERSION}'/g' -e 's/RELNUMBER/r'${REVISION}'/g' -e 's/DATE/'${DATE}'/g' > html/tfce.html

update: doc
	-@git fetch
	-@echo '% TFCE Toolbox' > Contents.m
	-@echo '% Version ' ${REVISION} ' ('${NEWVERSION}')' ${DATE} >> Contents.m
	-@cat Contents_info.txt >> Contents.m
	-@perl -p -i -e "s/${OLDVERSION}/${NEWVERSION}/g" spm_TFCE.m
	-@echo '% TFCE Toolbox' > INSTALL.txt
	-@echo '% Version ' ${REVISION} ${NEWVERSION} ${DATE} >> INSTALL.txt
	-@cat INSTALL_info.txt >> INSTALL.txt

zip: update
	-@echo zip
	-@test ! -d TFCE || rm -r TFCE
	-@mkdir TFCE
	-@cp -rp ${FILES} TFCE
	-@bash update_revision.sh
	-@zip ${ZIPFOLDER}/${ZIPFILE} -rm TFCE

scp: zip
	-@git tag -f ${VERSION} -m "Release ${VERSION}"
	-@echo scp to http://${STARGET_HOST}/tfce/${ZIPFILE}
	-@scp -O -P 2222 CHANGES.txt ${ZIPFOLDER}/${ZIPFILE} ${STARGET}
	-@bash -c "ssh -p 2222 ${STARGET_HOST} ln -fs ${STARGET_FOLDER}/${ZIPFILE} ${STARGET_FOLDER}/tfce_latest.zip"
	
