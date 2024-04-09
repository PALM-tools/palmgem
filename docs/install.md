# Installation
First, you need to install PostgreSQL database system with PostGIS gis extension.
To instal these systems, check and follow PostgreSQL and PostGIS documentation. On most Linux systems, these tools are part of the system distributions and can be install by standard Linux instalation tools.

## PostgreSQL + PostGIS
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


# Python libraries

How to install all python libraries and version. Different machines.
