#!/bin/bash
cd /home/oracle/scripts/kazutool
java ngsValidation/driver/DriverForGenerateGLstringFromSironaXML "/w_u/$1/upload/$2" "/tmp/$1.csv"
echo "SampleID	DB_Version	GL_String" > "/w_u/$1/upload/$3"
cat "/tmp/$1.csv" | sed -e 's/,/\t/g'  >> "/w_u/$1/upload/$3"
