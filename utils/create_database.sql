
\echo 'Connected to database'
select * from information_schema."information_schema_catalog_name";

\echo 'Register needed PALM-GEM functions';
\i palm_create_grid.sql;
\i palm_fill_building_holes.sql;
\i palm_fill_building_holes.sql;
commit;


