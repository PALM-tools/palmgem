#!/bin/bash

##############################
echo "This script MUST be run on PG server with user postgres or"
echo "another postgres user with privilege to create database."
echo "Press [Enter] key when ready to continue, press Ctrl+C to exit."
read anytext
echo -n "Enter database name: "; read dbname
echo -n "Enter user name or press [Enter] to use your name: "; read username
echo ""
##############################

if [ "x$dbname" == "x" ]; then
  echo "Use: create_database.sh <database_name> [<user_name>]"
  exit 1
fi

dbexists=$(psql -l | cut -d"|" -f1|grep -c "$dbname")
if [ ! $dbexists == 0 ]; then
  echo "Database $dbname just exists on the server!"
  echo "Use another database name or drop the existing database first."
  exit 1
fi


nr=$(psql -qtA -d postgres -c "select count(*) from pg_catalog.pg_roles where rolname = 'palmgem'")
if [ $nr == 0 ]; then
  echo "Create \"palmgem\" group role"
  psql -d postgres -c "create role \"palmgem\" createdb nologin"
fi

if [ "x$username" == "x" ]; then
  username=$USER
fi

nr=$(psql -qtA -d postgres -c "select count(*) from pg_catalog.pg_roles where rolname = '$username'")
if [ $nr == 0 ]; then
  echo "Create \"$username\" user role"
  echo -n "Enter password for user $username: "; read userpass
  psql -d postgres -c "create role \"$username\" login in role \"palmgem\""
  psql -d postgres -c "alter user \"$username\" with password '$userpass';"
fi

echo "Create database: $dbname, with owner \"palmgem\""
psql -d postgres -c "create database \"$dbname\" owner \"palmgem\";"
psql -d $dbname -c "create extension postgis;"
psql -d $dbname -c "create extension postgis_topology;"
psql -d $dbname -c "create extension intarray;"
psql -d $dbname -c "grant all on database \"$dbname\" to \"palmgem\" with grant option;"
psql -d $dbname -c "grant all on spatial_ref_sys to \"palmgem\";"
psql -d $dbname -c "grant all on spatial_ref_sys_srid_seq to \"palmgem\";"

echo "Database $dbname has been created and postgis has been enabled in it."
echo "Connect to the database as user $username and run SQL script create_database.sql."
echo "To run it you can use command:"
echo "psql -h <hostname> -p <port> -U <username> [-W] -d $dbname -f create_database.sql"
echo "Current connection info is:"
psql -d $dbname -c "\conninfo"
