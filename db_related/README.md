# IHIWS 17th Database
The IHIWS 17th database contains two main parts
1. the web application in https://ihiws17.stanford.edu/ 
2. the sftp server hidpl.stanford.edu

The functions are described in the paper and the scripts in this folder "db_related" could be used to re-build the database with minor changes.

## Web application
The web application was built with oracle SQL database (12c Standard Edition) and APEX 5.0. Before running the scripts, installation a compatible version of oracle SQL database and APEX is required.

After the installation and the successful login as a privileged user into Oracle APEX,

1. Edit the file oracle_sql/ihiws17_schema20171120.mssql so that the strings "/home/oracle/scripts/" are replaced with the folder containing scripts in the folder "scripts".
2. The web app allows users to execute pre-defined shell scripts with functions like merging multiple HML files, converting typing reports to IHIWS XML format with DBMS_SCHEDULER. Please set corresponding privileges for the schema holding the data.
3. Go to SQL Workshop -> SQL scripts -> Upload -> Choose File to upload oracle_sql/ihiws17_schema20171120.mssql and execute it to create the tables, functions, packages, views, procedures and etc used for the web app. 
4. Go to Application Builder -> Import -> Database Application, Page or Component Export to import the file oracle_sql/ihiws_apex_app.mssql and the application should be created.

## data

## scripts

## hlaPoly

Install package:
Change folder to db_related/
In R
>install.packages("hlaPoly",repos = NULL, type="source")



