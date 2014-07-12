#!/bin/bash

# create by lsw
# there are many problems of the arpack++ installed from the synaptic.
# this it the modified version created by simba.
# install arpack++ by running sudo ./install.sh
# see
# http://hi.baidu.com/europelee/item/38165ffbb52141683d148573

mkdir /usr/include/arpack++
cp ./include/* /usr/include/arpack++/
cd /usr/include/arpack++/
chmod a+r ./*