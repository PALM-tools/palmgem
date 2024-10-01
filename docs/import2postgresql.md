# Import dataset to PostgreSQL database
### Linux terminal version
Use prepared shell script in import_to_sql folder import_pgsql.sh to import downloaded data into postgresql database (please check names and paths of the importing files). After successful import you can check imported geodata using QGIS, detailed instruction can be found [here](visuallization.md).
### Python version (mainly for Windows)
Use prepared python script `import_sql.py`. Fill in a prepared configuration for importing `example_import.yaml` in `config` folder. Defined connection to database. Select if you want to create a schema from scratch (delete previous if exists, `scratch_import=True`) or append existing. Fill the path to files you would like to import. \
Run the script as follows:
```
python3 import_sql.py -c config_name.yaml
```
Log from the import is create based on import_schema name in `logs` folder.