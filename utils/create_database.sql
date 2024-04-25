/*
 * Copyright 2018-2024 Institute of Computer Science of the Czech Academy of
 * Sciences, Prague, Czech Republic. Authors: Martin Bures, Jaroslav Resler.
 *
 * This file is part of PALM-GeM.
 *
 * PALM-GeM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * PALM-GeM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * PALM-GeM. If not, see <https://www.gnu.org/licenses/>.
 */

\echo 'Connected to database'
select * from information_schema."information_schema_catalog_name";

\echo 'Register needed PALM-GEM functions';
\i palm_create_grid.sql;
\i palm_fill_building_holes.sql;
\i palm_fill_building_holes.sql;
commit;


