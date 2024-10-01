# Installation
First, you need to install PostgreSQL database system with PostGIS gis extension.
To instal these systems, check and follow PostgreSQL and PostGIS documentation. On most Linux systems, these tools are part of the system distributions and can be install by standard Linux instalation tools.

## PostgreSQL + PostGIS, linux system
Brief configuration steps after the instalation of PostgreSQL and PostGIS:

If not already started, start (possibly as root) the postgres sql server.

_systemctl start postgresql_

Create database:

  - You have to either copy intialization scripts in directory
  _utils_ somewhere where user postgres can read it, e.g.

_cp -ra utils /var/lib/pgsql/_

  - Switch to user 'postgres:

_sudo -u postgres -i_

or if you don't use sudo:

-su; su - postgres_

  - run the script

_util/create_database.sh_

Alternatively instead of running the script, you can follow individual tasks in it step by step.

Logout from the database, the rest of the configuration is done under your account.

Loading required SQL functions into your database:

_cd utils_

_psql -h hostname -p port -U username [-W] -d dbname -f create_database.sql_

## PostgreSQL + PostGIS, Windows system
1) Download postgresql database from their website https://www.enterprisedb.com/downloads/postgres-postgresql-downloads
2) Install postgreSQL from downloaded .exe file
3) Setup password for access
4) Run pgadmin4 
   1) Create a group role: Login/Group role
   2) Name it palmgem
   3) Add all privileges
5) Create a user role: Login/Group role 
   1) Name it according to user
   2) Add sufficient privileges (if you are only user, add all)
6) Create Database with name, e.g. palm_static
   1) Owner to palmgem
7) Install Postgis extension
   1) Use Stack Builder and tutorial here: https://www.bostongis.com/PrinterFriendly.aspx?content_name=postgis_tut01
   2) Keep information about installation location
   3) In query tool in palm_static database execute:
      1) `create extension postgis`
      2) `create extension postgis_topology`
      3) `create extension intarray`
      4) `create extension postgis_raster`
   4) Change owner ship of palm_static database to palmgem superuser.
   5) Execute in query tool in palm_static database
      1) `grant all on database palm_static to palmgem with grant option`
      2) `grant all on spatial_ref_sys to palmgem`
      3) (dont, need testing) `grant all on spatial_ref_sys_srid_seq to palmgem `
8) If you share the database with colleagues, create a common group with permissions for all colleagues, e.g. named palm
9) Upload all sql function `palm_create_grid.sql`, `palm_fill_building_holes.sql`, `palm_surface.sql` and `palm_tree_grid.sql` by copying text from function into query builder and execute

## Python libraries
PALM-GeM may be used out-of-the box with the project directory as long as the all
required libraries are available. The easiest way to install them is using

```
pip3 install -r requirements.txt
```

however if you prefer slightly different versions of the libraries specified in
`requirements.txt`, potentially from your operating system's distribution, you
may try them as well.