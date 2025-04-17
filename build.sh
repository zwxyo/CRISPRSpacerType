#!/bin/bash


mkdir -p $PREFIX/share/CRISPRSpacerType
cp -r src config database bin $PREFIX/share/CRISPRSpacerType

mkdir -p $PREFIX/bin
chmod +x $PREFIX/share/CRISPRSpacerType/bin/*.sh

ln -s $PREFIX/share/CRISPRSpacerType/bin/CRISPRSpacerType.sh $PREFIX/bin/CRISPRSpacerType
