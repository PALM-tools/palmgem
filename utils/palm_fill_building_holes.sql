CREATE SEQUENCE IF NOT EXISTS spatial_ref_sys_srid_seq MINVALUE 100000;

/************************
* Function that fills holes in building. Holes in building are defined that >= 3 neighbors are defined as building while
* hole grid is defined as non building
*************************/
drop function if exists palm_fill_building_holes(
    case_schema text,
    grid_table text,
    landcover_table text,
    building_min integer,
    building_max integer,
    nx integer,
    ny integer,
    debug_level integer);

create or replace function palm_fill_building_holes(
    case_schema text,
    grid_table text,
    landcover_table text,
    building_min integer,
    building_max integer,
    nx integer,
    ny integer,
    debug_level integer)
  returns boolean as
$$
declare
    ret boolean;
    res text;
    sqltext text;
	sqlcount text;
    sqlfill text;
    i integer;
    j integer;
	count_fill integer;
	iters integer;
	max_count integer;
begin
    ret = false;
    sqltext = format('ALTER TABLE %I.%I DROP CONSTRAINT IF EXISTS landcover_pkey',case_schema, landcover_table);
    execute(sqltext);
    sqltext = format('ALTER TABLE %I.%I ADD PRIMARY KEY (lid)',case_schema, landcover_table);
    execute(sqltext);
    sqltext = format('ALTER TABLE %I.%I DROP CONSTRAINT IF EXISTS grid_pkey',case_schema, grid_table);
    execute(sqltext);
    sqltext = format('ALTER TABLE %I.%I ADD PRIMARY KEY (i,j,lid)',case_schema, grid_table);
    execute(sqltext);

	sqlcount = format('SELECT COUNT(*) FROM %I.%I AS g ' ||
                    'LEFT OUTER JOIN %I.%I AS l ON g.lid=l.lid ' ||
                    'WHERE ((SELECT COUNT(*) FROM %I.%I AS b ' ||
	   				'	        LEFT OUTER JOIN %I.%I AS ll ON b.lid=ll.lid ' ||
	   				' 		         WHERE  ((b.i = g.i+1 and b.j = g.j) OR ' ||
                    '                       (b.i = g.i-1 and b.j = g.j) OR ' ||
                    '                       (b.i = g.i   and b.j = g.j-1) OR ' ||
                    '                       (b.i = g.i   and b.j = g.j+1)) and ll.type BETWEEN $1 AND $2 ' ||
                    '       ) + CASE WHEN ((g.i = 0 OR g.i = $3) AND l.type NOT BETWEEN $1 AND $2) THEN 1 ELSE 0 END ' ||
                    '         + CASE WHEN ((g.j = 0 OR g.j = $4) AND l.type NOT BETWEEN $1 AND $2) THEN 1 ELSE 0 END '
                    '       ) > 2 AND l.type NOT BETWEEN $1 AND $2 ', case_schema, grid_table,
                                                                        case_schema, landcover_table,
                                                                        case_schema, grid_table,
                                                                        case_schema, landcover_table);
--    sqltext = format('SELECT g.i, g.j FROM %I.%I AS g ' ||
--                    'LEFT OUTER JOIN %I.%I AS l ON g.lid=l.lid ' ||
--                    'WHERE (SELECT COUNT(*) FROM %I.%I AS b ' ||
--	   				'	        LEFT OUTER JOIN %I.%I AS ll ON b.lid=ll.lid ' ||
--	   				' 		         WHERE  ((b.i = g.i+1 and b.j = g.j) OR ' ||
--                    '                       (b.i = g.i-1 and b.j = g.j) OR ' ||
--                    '                       (b.i = g.i   and b.j = g.j-1) OR ' ||
--                    '                       (b.i = g.i   and b.j = g.j+1)) and ll.type BETWEEN 900 AND 999 ' ||
--                    '       ) > 2 AND l.type NOT BETWEEN 900 AND 999 ', case_schema, grid_table,
--                                                                        case_schema, landcover_table,
--                                                                        case_schema, grid_table,
--                                                                        case_schema, landcover_table);

    sqltext = format('SELECT g.i, g.j FROM %I.%I AS g ' ||
                    'LEFT OUTER JOIN %I.%I AS l ON g.lid=l.lid ' ||
                    'WHERE ((SELECT COUNT(*) FROM %I.%I AS b ' ||
	   				'	        LEFT OUTER JOIN %I.%I AS ll ON b.lid=ll.lid ' ||
	   				' 		         WHERE  ((b.i = g.i+1 and b.j = g.j) OR ' ||
                    '                       (b.i = g.i-1 and b.j = g.j) OR ' ||
                    '                       (b.i = g.i   and b.j = g.j-1) OR ' ||
                    '                       (b.i = g.i   and b.j = g.j+1)) and ll.type BETWEEN $1 AND $2 ' ||
                    '       ) + CASE WHEN ((g.i = 0 OR g.i = $3) AND l.type NOT BETWEEN $1 AND $2) THEN 1 ELSE 0 END ' ||
                    '         + CASE WHEN ((g.j = 0 OR g.j = $4) AND l.type NOT BETWEEN $1 AND $2) THEN 1 ELSE 0 END '
                    '       ) > 2 ' ||
                    ' AND l.type NOT BETWEEN $1 AND $2 ', case_schema, grid_table,
                                                                        case_schema, landcover_table,
                                                                        case_schema, grid_table,
                                                                        case_schema, landcover_table);

    sqlfill = format('UPDATE %I.%I AS g SET lid = (SELECT b.lid FROM %I.%I AS b ' ||
                        'LEFT OUTER JOIN %I.%I AS l ON b.lid=l.lid ' ||
                            'WHERE  (((b.i = g.i+1 and b.j = g.j) OR ' ||
                    '                (b.i = g.i-1 and b.j = g.j) OR ' ||
                    '                (b.i = g.i   and b.j = g.j-1) OR ' ||
                    '                (b.i = g.i   and b.j = g.j+1)) and l.type BETWEEN $1 AND $2) ' ||
								   'GROUP BY b.lid ORDER BY COUNT(*) LIMIT 1) ' ||
	                    'WHERE g.i = $3 AND g.j = $4', case_schema, grid_table,
	                    case_schema, grid_table, case_schema, landcover_table);

	execute(sqlcount) using building_min, building_max, (nx-1), (ny-1) into count_fill;
	if debug_level < 4 then
	    raise notice 'In 0 iteration, % grids need to be filled', count_fill;
	end if;

	iters = 0;
	--TODO: introduce max loop, just in case of infinite cycle
	max_count = 20;
	while count_fill > 0 loop
		iters = iters + 1;
		for i, j in execute sqltext using building_min, building_max, (nx-1), (ny-1) loop
			execute sqlfill using building_min, building_max, i, j;
			--raise notice 'grid [i,j]=[%,%] has been modified to building', i, j;
		end loop;
		execute(sqlcount) using building_min, building_max, (nx-1), (ny-1) into count_fill;
		if debug_level < 3 then
	        raise notice '%   In % iteration, % grids need to be filled', clock_timestamp()::timestamp(0), iters, count_fill;
	    end if;
	    if iters > max_count then
	        raise notice 'Max count has been achieved, stop the loop';
	    end if;

	end loop;

    --sqltext = format('ALTER TABLE %I.%I DROP CONSTRAINT grid_pkey',case_schema, grid_table);
    --execute(sqltext);
    --sqltext = format('ALTER TABLE %I.%I ADD PRIMARY KEY (id)',case_schema, grid_table);
    --execute(sqltext);

    ret = true;
    return ret;
end
$$
language plpgsql volatile
cost 100;