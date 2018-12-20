# Personal Makefile variables
#
# $Id$

VERSION=`svn info |grep Revision|sed -e 's/Revision: //g'`
DATE=`svn info |grep 'Last Changed Date: '|sed -e 's/Last Changed Date: //g'|cut -f1 -d' '`

TARGET=/Users/gaser/spm/spm12/toolbox/TFCE
TARGET2=/Volumes/UltraMax/spm12/toolbox/TFCE

STARGET_HOST=dbm.neuro.uni-jena.de
STARGET_HTDOCS=${STARGET_HOST}:/Applications/xampp/htdocs/
STARGET_FOLDER=/Applications/xampp/htdocs/tfce
STARGET=${STARGET_HOST}:${STARGET_FOLDER}

MATLAB_FILES=Contents.m cg_get_tfce_results.m cg_tfce_results.m cg_tfce_list.m cg_tfce_update.m cg_progress.m spm_TFCE.m tfce_mesh.m snpm_P_FDR.m tbx_cfg_tfce_estimate.m cg_tfce_estimate.m cg_tfce_surf_max.m
C_FILES=tfceMex_pthread.* 
MISC_FILES=TFCE.man

FILES=${MATLAB_FILES} ${C_FILES} ${MISC_FILES}

ZIPFILE=tfce_r$(VERSION).zip

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

help:
	-@echo Available commands:
	-@echo install zip scp update

update:
	-@svn update
	-@echo '% TFCE Toolbox' > Contents.m
	-@echo '% Version ' ${VERSION} ' (TFCE) ' ${DATE} >> Contents.m
	-@cat Contents_info.txt >> Contents.m
	-@echo '% __________________________________________________________________________' > TFCE.man
	-@echo '% TFCE Toolbox for SPM8/SPM12' >> TFCE.man
	-@echo '% Version ' ${VERSION} ' (TFCE) ' ${DATE} >> TFCE.man
	-@cat TFCE.txt >> TFCE.man
	-@echo '% TFCE Toolbox' > INSTALL.txt
	-@echo '% Version ' ${VERSION} ' (TFCE) ' ${DATE} >> INSTALL.txt
	-@cat INSTALL_info.txt >> INSTALL.txt

zip: update
	-@echo zip
	-@test ! -d TFCE || rm -r TFCE
	-@mkdir TFCE
	-@cp -rp ${FILES} TFCE
	-@zip ${ZIPFILE} -rm TFCE

scp: zip
	-@echo scp to http://dbm.neuro.uni-jena.de/tfce/${ZIPFILE}
	-@scp -P 2222 CHANGES.txt ${ZIPFILE} ${STARGET}
	-@bash -c "ssh ${STARGET_HOST} ln -Fs ${STARGET_FOLDER}/${ZIPFILE} ${STARGET_FOLDER}/tfce_latest.zip"
