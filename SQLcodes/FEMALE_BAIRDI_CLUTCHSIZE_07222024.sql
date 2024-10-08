-- This script produces a table of population estimates for bairdi Tanner
-- crab from the 1975-2016 EBS trawl surveys.  Population is calculated
-- for males and female crab by stock assessment size category or maturity for each district. 
-- Abundance is calculated for unsexed crab, but biomass is not as we have no
-- size-weight regression factors for unsexed crab.

-- This script requires as input the master crab table (ebscrab) populated
-- with the survey data to be analyzed; a subset of the racebase.haul
-- table containing the haul data for the cruises being analyzed; 
-- and a strata lookup table.
-------------------------------------------------------------------------------

-- Don't want to average non-rkc in BB retow years, get rid of haul type 17 tows

drop table haul_newtimeseries_noretow;

create table haul_newtimeseries_noretow as
select * from haul_newtimeseries
where haul_type <> 17;


-- Create tables of raw catch by 1-mm size bin and sex
-- Separate by sex because male size group categories require shell condition
-- and female weights (post-2009) require clutch size



-- Females (done separately from males because need clutch size info)

drop table cb_number_size1_female;

create table cb_number_size1_female as
select c.hauljoin,c.vessel,c.cruise,c.haul,h.gis_station,species_code,shell_condition, clutch_size,
(trunc(width/1) * 1)size1,
(sum(CASE
		 when species_code = 68560
		 and sex = 2
		 then sampling_factor
		 else 0
		 end)) number_female_size1
from crab.ebscrab c, haul_newtimeseries_noretow h
where species_code = 68560
and width <> 999
and c.hauljoin(+) = h.hauljoin
and haul_type <> 17
group by c.hauljoin,
	  	 c.vessel,
		 c.cruise,
		 c.haul,
		 h.gis_station,
		 species_code,
         shell_condition,
         clutch_size,
		 (trunc(width/1) * 1);


-- unsexed



--  This section calculates the weight of the bairdi Tanner crab by haul, sex,
--  shell condition and 1-mm size group.  A width-weight regression
--  factor is applied, and multiplied by the number of crab caught in that
--  haul/sex/shellcon/size bin (from above section).  
--  The regression factor does not include unsexed crab, therefore no weights
--  will be calculated for unsexed crab




drop table cb_weight_grams_female;

create table cb_weight_grams_female as
select hauljoin,vessel,cruise,haul,gis_station,species_code,shell_condition, clutch_size,size1,
(CASE
--    WHEN cruise < 201001
--      THEN ((0.00182 * (power(size1,2.70462))) * number_female_size1)
    WHEN cruise >= 197501 and clutch_size <= 1
      THEN ((0.000562 * (power(size1,2.816928))) * number_female_size1)
    WHEN cruise >= 197501 and clutch_size > 1
      THEN ((0.000441 * (power(size1,2.898686))) * number_female_size1)
    ELSE 0
    END) wgt_female_size1
from cb_number_size1_female
order by cruise,vessel,haul,gis_station,size1;

-- Using actual female maturity in this run, so select for clutch size here

drop table cb_number_size1_matfem;

create table cb_number_size1_matfem as
select hauljoin, vessel, cruise, haul, gis_station, species_code,shell_condition, clutch_size, size1,
	   (sum(CASE
	   			WHEN   clutch_size = 0  
				   THEN number_female_size1
				ELSE 0
				END))  number_female_size1_immature,	
	   (sum(CASE
	   			WHEN   clutch_size >= 1  
                   and clutch_size <999
				   THEN number_female_size1
				ELSE 0
				END))  number_female_size1_mature
    from cb_number_size1_female
    where species_code = 68560
		 --and width <> 999
	     --and hauljoin(+) = h.hauljoin
       group by hauljoin,
	            vessel,
				cruise,
				haul,
				gis_station,
				species_code,
                shell_condition,
                clutch_size,
				size1;
        
-- And calculate weight of crab by actual maturity        

drop table cb_weight_grams_matfem;

create table cb_weight_grams_matfem as
select hauljoin, vessel, cruise, haul, gis_station, species_code, shell_condition,clutch_size, size1,
	   (sum(CASE
	   			WHEN   clutch_size = 0  
				   THEN wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_immature,	
	   (sum(CASE
	   			WHEN   clutch_size >= 1
                   and clutch_size <999
				   THEN wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_mature
    from cb_weight_grams_female
    where species_code = 68560
		 --and width <> 999
	     --and hauljoin(+) = h.hauljoin
       group by hauljoin,
	            vessel,
				cruise,
				haul,
				gis_station,
				species_code,
                shell_condition,
                clutch_size,
				size1;

-- combine male and female weight tables

drop table cb_weight_grams_size1;

create table cb_weight_grams_size1 
( HAULJOIN                        NUMBER(12),
  VESSEL                          NUMBER(4),
  CRUISE                          NUMBER(6),
  HAUL                            NUMBER(4),
  GIS_STATION                     VARCHAR2(10),
  SPECIES_CODE                    NUMBER(6),
  SHELL_CONDITION                 NUMBER,
  CLUTCH_SIZE                     NUMBER,
  SIZE1                           NUMBER,
  WGT_FEMALE_SIZE1_IMMATURE       NUMBER,
  WGT_FEMALE_SIZE1_MATURE         NUMBER
);


insert into cb_weight_grams_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,shell_condition,clutch_size,size1,
wgt_female_size1_immature,wgt_female_size1_mature
from cb_weight_grams_matfem;

-- convert to metric tons

drop table cb_weight_mt_size1;

create table cb_weight_mt_size1 as
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,shell_condition,clutch_size,
(wgt_female_size1_immature * 0.000001) mt_female_size1_immature,
(wgt_female_size1_mature * 0.000001) mt_female_size1_mature
from cb_weight_grams_size1
order by cruise,vessel,haul,gis_station,size1;

-- Combine the male, female, and unsexed by number tables

drop table cb_number_size1;

create table cb_number_size1 
( HAULJOIN                        NUMBER(12),
  VESSEL                          NUMBER(4),
  CRUISE                          NUMBER(6),
  HAUL                            NUMBER(4),
  GIS_STATION                     VARCHAR2(10),
  SPECIES_CODE                    NUMBER(6),
  SHELL_CONDITION                 NUMBER,
  CLUTCH_SIZE                     NUMBER,
  SIZE1                           NUMBER,
  NUMBER_FEMALE_SIZE1_IMMATURE    NUMBER,
  NUMBER_FEMALE_SIZE1_MATURE      NUMBER
);


insert into cb_number_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,shell_condition, clutch_size, size1,
number_female_size1_immature,number_female_size1_mature
from cb_number_size1_matfem;


-- This section sums the bairdi Tanner crab catch records by haul, sex,
-- and 1-mm size group.  

drop table cb_number_sizegroup;

create table cb_number_sizegroup as
select hauljoin, vessel, cruise, haul, gis_station, species_code, 
	   
	   (sum(CASE
	   			WHEN  clutch_size = 1  
			    THEN number_female_size1_mature
				ELSE 0
				END))  number_female_barren,
     (sum(CASE
	   			WHEN clutch_size = 2
			    THEN number_female_size1_mature
				ELSE 0
				END))  number_female_trace,	                  
                
	   (sum(CASE
	   			WHEN  clutch_size = 3
			    THEN number_female_size1_mature
				ELSE 0
				END))  number_female_quarter,	
				
	   (sum(CASE
	   			WHEN  clutch_size = 4
			    THEN number_female_size1_mature
				ELSE 0
				END))  number_female_half,	          
	   
	   (sum(CASE
	   			WHEN  clutch_size = 5
			    THEN number_female_size1_mature
				ELSE 0
				END))  number_female_three_quarter,
        
	   (sum(CASE
	   			WHEN  clutch_size = 6
			    THEN number_female_size1_mature
				ELSE 0
				END))  number_female_full,	        	
	   sum(number_female_size1_immature) number_female_immature,
     sum(number_female_size1_mature) number_female_mature,
     (sum(number_female_size1_immature)+ sum(number_female_size1_mature)) number_female_total												
    from cb_number_size1
    where species_code = 68560
    group by hauljoin,
	            vessel,
				cruise,
				haul,
				gis_station,
				species_code;
				

drop table cb_weight_mt_sizegroup;

create table cb_weight_mt_sizegroup as
select hauljoin, vessel, cruise, haul, gis_station, species_code, 
	   (sum(CASE
	   			WHEN  clutch_size = 1
			    THEN mt_female_size1_mature
				ELSE 0
				END))  mt_female_barren,
      (sum(CASE
	   			WHEN  clutch_size = 2
			    THEN mt_female_size1_mature
				ELSE 0
				END))  mt_female_trace,
      (sum(CASE
	   			WHEN  clutch_size = 3
			    THEN mt_female_size1_mature
				ELSE 0
				END))  mt_female_quarter,
      (sum(CASE
	   			WHEN  clutch_size = 4
			    THEN mt_female_size1_mature
				ELSE 0
				END))  mt_female_half,
      (sum(CASE
	   			WHEN  clutch_size = 5
			    THEN mt_female_size1_mature
				ELSE 0
				END))  mt_female_three_quarter,
      (sum(CASE
	   			WHEN  clutch_size = 6
			    THEN mt_female_size1_mature
				ELSE 0
				END))  mt_female_full,
     sum(mt_female_size1_immature) mt_female_immature,
     sum(mt_female_size1_mature) mt_female_mature,
     (sum(mt_female_size1_immature)+ sum(mt_female_size1_mature)) mt_female_total
	   from cb_weight_mt_size1
         where species_code = 68560
       group by hauljoin,
	            vessel,
				cruise,
				haul,
				gis_station,
				species_code;
				
				

-- This section combines the haul and catch data, including
-- those haul/size groups where there was no catch.				

drop table cb_num_sizegroup_union;

create table cb_num_sizegroup_union as
select h.hauljoin,h.vessel,h.cruise,h.haul,mid_latitude,mid_longitude,h.gis_station,survey_year,
nvl(species_code,68560) species_code,
nvl(number_female_barren,0) number_female_barren,
nvl(number_female_trace,0) number_female_trace,
nvl(number_female_quarter,0) number_female_quarter,
nvl(number_female_half,0) number_female_half,
nvl(number_female_three_quarter,0) number_female_three_quarter,
nvl(number_female_full,0) number_female_full,
nvl(number_female_immature,0) number_female_immature,
nvl(number_female_mature,0) number_female_mature,
nvl(number_female_total,0) number_female_total
from haul_newtimeseries_noretow h full outer join cb_number_sizegroup c
on h.hauljoin = c.hauljoin
where haul_type <> 17;

--  Similarly, by weight.

drop table cb_wgt_sizegroup_union;

create table cb_wgt_sizegroup_union as
select h.hauljoin,h.vessel,h.cruise,h.haul,mid_latitude,mid_longitude,h.gis_station,survey_year,
nvl(species_code,68560) species_code,
nvl(mt_female_barren,0) mt_female_barren,
nvl(mt_female_trace,0) mt_female_trace,
nvl(mt_female_quarter,0) mt_female_quarter,
nvl(mt_female_half,0) mt_female_half,
nvl(mt_female_three_quarter,0) mt_female_three_quarter,
nvl(mt_female_full,0) mt_female_full,
nvl(mt_female_immature,0) mt_female_immature,
nvl(mt_female_mature,0) mt_female_mature,
nvl(mt_female_total,0) mt_female_total
from haul_newtimeseries_noretow h full outer join cb_weight_mt_sizegroup c
on h.hauljoin = c.hauljoin
where haul_type <> 17;


-- This section calculates cpue for each haul.
-- If a station contains multiple tows, cpue
-- is calculated for each of the tows, not averaged for the station.
-- A value, even if 0 for no catch, is output for every size group,
-- every haul.  CPUE is calculated as number of crabs per square
-- nautical mile towed; area swept is the distance fished multiplied
-- by the actual (measured) net width.

drop table cb_cpuenum_sizegroup;

create table cb_cpuenum_sizegroup as
select c.hauljoin,c.vessel,c.cruise,c.haul,h.mid_latitude,h.mid_longitude,c.gis_station,c.survey_year,
c.species_code,(((h.net_width/1000) * h.distance_fished) * 0.29155335) area_swept,
(number_female_barren / (((net_width/1000) * distance_fished) * 0.29155335)) female_barren_cpuenum,
(number_female_trace / (((net_width/1000) * distance_fished) * 0.29155335)) female_trace_cpuenum,
(number_female_quarter / (((net_width/1000) * distance_fished) * 0.29155335)) female_quarter_cpuenum,
(number_female_half / (((net_width/1000) * distance_fished) * 0.29155335)) female_half_cpuenum,
(number_female_three_quarter / (((net_width/1000) * distance_fished) * 0.29155335)) female_three_quarter_cpuenum,
(number_female_full / (((net_width/1000) * distance_fished) * 0.29155335)) female_full_cpuenum,
(number_female_immature / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_immature,
(number_female_mature / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_mature,
(number_female_total / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_total
from cb_num_sizegroup_union c, haul_newtimeseries_noretow h
where c.hauljoin = h.hauljoin
and haul_type <> 17;


-- This section calculates cpue by weight for each haul.
-- If a station contains multiple tows, cpue is calculated 
-- for each of the tows, not averaged for the station.
-- A value, even if 0 for no catch, is output for every size group,
-- every haul.  CPUE is calculated as weight of crabs (already converted to metric tons) 
-- per square nautical mile towed; area swept is the distance fished multiplied
-- by the actual (measured) net width.

drop table cb_cpuewgt_sizegroup;

create table cb_cpuewgt_sizegroup as
select c.hauljoin,c.vessel,c.cruise,c.haul,h.mid_latitude,h.mid_longitude,
c.gis_station,c.survey_year,c.species_code,(((h.net_width/1000) * h.distance_fished) * 0.29155335) area_swept,
(mt_female_barren / (((net_width/1000) * distance_fished) * 0.29155335)) female_barren_cpuewgt,
(mt_female_trace / (((net_width/1000) * distance_fished) * 0.29155335)) female_trace_cpuewgt,
(mt_female_quarter / (((net_width/1000) * distance_fished) * 0.29155335)) female_quarter_cpuewgt,
(mt_female_half / (((net_width/1000) * distance_fished) * 0.29155335)) female_half_cpuewgt,
(mt_female_three_quarter / (((net_width/1000) * distance_fished) * 0.29155335)) female_three_quarter_cpuewgt,
(mt_female_full / (((net_width/1000) * distance_fished) * 0.29155335)) female_full_cpuewgt,
(mt_female_immature / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_immature,
(mt_female_mature / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_mature,
(mt_female_total / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_total
from cb_wgt_sizegroup_union c, haul_newtimeseries_noretow h
where c.hauljoin = h.hauljoin
and haul_type <> 17;


drop table cb_meancpuenum_sizegroup;

create table cb_meancpuenum_sizegroup as
select c.survey_year,district,
AVG (female_barren_cpuenum) meancpuenum_female_barren,
AVG (female_trace_cpuenum) meancpuenum_female_trace,
AVG (female_quarter_cpuenum) meancpuenum_female_quarter,
AVG (female_half_cpuenum) meancpuenum_female_half,
AVG (female_three_quarter_cpuenum) meancpuenum_female_three_quarter,
AVG (female_full_cpuenum) meancpuenum_female_full,
AVG (female_cpuenum_immature) meancpuenum_female_immature,
AVG (female_cpuenum_mature) meancpuenum_female_mature,
AVG (female_cpuenum_total) meancpuenum_female_total
from cb_cpuenum_sizegroup c, strata_bairdi_newtimeseries s
where c.gis_station = s.station_id
and c.survey_year = s.survey_year
group by c.survey_year,district;

drop table cb_meancpuewgt_sizegroup;

create table cb_meancpuewgt_sizegroup as
select c.survey_year,district,
AVG (female_barren_cpuewgt) meancpuewgt_female_barren,
AVG (female_trace_cpuewgt) meancpuewgt_female_trace,
AVG (female_quarter_cpuewgt) meancpuewgt_female_quarter,
AVG (female_half_cpuewgt) meancpuewgt_female_half,
AVG (female_three_quarter_cpuewgt) meancpuewgt_female_three_quarter,
AVG (female_full_cpuewgt) meancpuewgt_female_full,
AVG (female_cpuewgt_immature) meancpuewgt_female_immature,
AVG (female_cpuewgt_mature) meancpuewgt_female_mature,
AVG (female_cpuewgt_total) meancpuewgt_female_total
from cb_cpuewgt_sizegroup c, strata_bairdi_newtimeseries s
where c.gis_station = s.station_id
and c.survey_year = s.survey_year
group by c.survey_year,district;


drop table cb_popbystratum_sizegroup;

create table cb_popbystratum_sizegroup as
select distinct c.survey_year,stratum,c.district,
 (meancpuenum_female_barren * total_area) pop_female_barren,
 (meancpuenum_female_trace * total_area) pop_female_trace,
 (meancpuenum_female_quarter * total_area) pop_female_quarter,
 (meancpuenum_female_half * total_area) pop_female_half,
 (meancpuenum_female_three_quarter * total_area) pop_female_three_quarter,
 (meancpuenum_female_full * total_area) pop_female_full,
 (meancpuenum_female_immature * total_area) pop_female_immature,
 (meancpuenum_female_mature * total_area) pop_female_mature,
 (meancpuenum_female_total * total_area) pop_female_total
from cb_meancpuenum_sizegroup c, strata_bairdi_newtimeseries s
where c.district = s.district
and c.survey_year = s.survey_year
order by survey_year,district;

drop table cb_biobystratum_sizegroup;

create table cb_biobystratum_sizegroup as
select distinct c.survey_year,stratum,c.district,
(meancpuewgt_female_barren * total_area) bio_female_barren,
(meancpuewgt_female_trace * total_area) bio_female_trace,
(meancpuewgt_female_quarter * total_area) bio_female_quarter,
(meancpuewgt_female_half * total_area) bio_female_half,
(meancpuewgt_female_three_quarter * total_area) bio_female_three_quarter,
(meancpuewgt_female_full * total_area) bio_female_full,
(meancpuewgt_female_immature * total_area) bio_female_immature,
(meancpuewgt_female_mature * total_area) bio_female_mature,
(meancpuewgt_female_total * total_area) bio_female_total
from cb_meancpuewgt_sizegroup c, strata_bairdi_newtimeseries s
where c.district = s.district
and c.survey_year = s.survey_year
order by survey_year,district;


drop table cb_varcpuenum_sizegroup;

create table cb_varcpuenum_sizegroup as
select c.survey_year,district,
VARIANCE (female_barren_cpuenum) varcpuenum_female_barren,
VARIANCE (female_trace_cpuenum) varcpuenum_female_trace,
VARIANCE (female_quarter_cpuenum) varcpuenum_female_quarter,
VARIANCE (female_half_cpuenum) varcpuenum_female_half,
VARIANCE (female_three_quarter_cpuenum) varcpuenum_female_three_quarter,
VARIANCE (female_full_cpuenum) varcpuenum_female_full,
VARIANCE (female_cpuenum_immature) varcpuenum_female_immature,
VARIANCE (female_cpuenum_mature) varcpuenum_female_mature,
VARIANCE (female_cpuenum_total) varcpuenum_female_total
from cb_cpuenum_sizegroup c, strata_bairdi_newtimeseries s
where c.gis_station = s.station_id
and c.survey_year = s.survey_year
group by c.survey_year,district;

drop table cb_varcpuewgt_sizegroup;

create table cb_varcpuewgt_sizegroup as
select c.survey_year,district,
VARIANCE (female_barren_cpuewgt) varcpuewgt_female_barren,
VARIANCE (female_trace_cpuewgt) varcpuewgt_female_trace,
VARIANCE (female_quarter_cpuewgt) varcpuewgt_female_quarter,
VARIANCE (female_half_cpuewgt) varcpuewgt_female_half,
VARIANCE (female_three_quarter_cpuewgt) varcpuewgt_female_three_quarter,
VARIANCE (female_full_cpuewgt) varcpuewgt_female_full,
VARIANCE (female_cpuewgt_immature) varcpuewgt_female_immature,
VARIANCE (female_cpuewgt_mature) varcpuewgt_female_mature,
VARIANCE (female_cpuewgt_total) varcpuewgt_female_total
from cb_cpuewgt_sizegroup c, strata_bairdi_newtimeseries s
where c.gis_station = s.station_id
and c.survey_year = s.survey_year
group by c.survey_year,district;


drop table cb_haulcount;

create table cb_haulcount as
select count(hauljoin)number_tows, h.survey_year, district
from haul_newtimeseries_noretow h, strata_bairdi_newtimeseries s
where h.gis_station = s.station_id
and h.survey_year = s.survey_year
and haul_type <> 17
group by h.survey_year, district;

drop table cb_variancepop_sizegroup;

create table cb_variancepop_sizegroup as
select distinct c.survey_year,stratum,c.district,
((varcpuenum_female_barren * (power(total_area,2)))/number_tows) varpop_female_barren,
((varcpuenum_female_trace * (power(total_area,2)))/number_tows) varpop_female_trace,
((varcpuenum_female_quarter * (power(total_area,2)))/number_tows) varpop_female_quarter,
((varcpuenum_female_half * (power(total_area,2)))/number_tows) varpop_female_half,
((varcpuenum_female_three_quarter * (power(total_area,2)))/number_tows) varpop_female_three_quarter,
((varcpuenum_female_full * (power(total_area,2)))/number_tows) varpop_female_full,
((varcpuenum_female_immature * (power(total_area,2)))/number_tows) varpop_female_immature,
((varcpuenum_female_mature * (power(total_area,2)))/number_tows) varpop_female_mature,
((varcpuenum_female_total * (power(total_area,2)))/number_tows) varpop_female_total
from strata_bairdi_newtimeseries s, cb_varcpuenum_sizegroup c, cb_haulcount n
where c.district = s.district
and c.district = n.district
and c.survey_year = s.survey_year
and c.survey_year = n.survey_year
order by c.survey_year,stratum;

drop table cb_variancebio_sizegroup;

create table cb_variancebio_sizegroup as
select distinct c.survey_year,stratum,c.district,
((varcpuewgt_female_barren * (power(total_area,2)))/number_tows) varbio_female_barren,
((varcpuewgt_female_trace * (power(total_area,2)))/number_tows) varbio_female_trace,
((varcpuewgt_female_quarter * (power(total_area,2)))/number_tows) varbio_female_quarter,
((varcpuewgt_female_half * (power(total_area,2)))/number_tows) varbio_female_half,
((varcpuewgt_female_three_quarter * (power(total_area,2)))/number_tows) varbio_female_three_quarter,
((varcpuewgt_female_full * (power(total_area,2)))/number_tows) varbio_female_full,
((varcpuewgt_female_immature * (power(total_area,2)))/number_tows) varbio_female_immature,
((varcpuewgt_female_mature * (power(total_area,2)))/number_tows) varbio_female_mature,
((varcpuewgt_female_total * (power(total_area,2)))/number_tows) varbio_female_total
from strata_bairdi_newtimeseries s, cb_varcpuewgt_sizegroup c, cb_haulcount n
where c.district = s.district
and c.district = n.district
and c.survey_year = s.survey_year
and c.survey_year = n.survey_year
order by c.survey_year,stratum;


-- Calculation by stock or district from this point on
-- For bairdi, this will be total, east of 166W, and west of 166W

drop table cb_popall_sizegroup;

create table cb_popall_sizegroup as
select survey_year,
sum(pop_female_barren) sum_pop_female_barren,
sum(pop_female_trace) sum_pop_female_trace,
sum(pop_female_quarter) sum_pop_female_quarter,
sum(pop_female_half) sum_pop_female_half,
sum(pop_female_three_quarter) sum_pop_female_three_quarter,
sum(pop_female_full) sum_pop_female_full,
sum(pop_female_immature) sum_pop_female_immature,
sum(pop_female_mature) sum_pop_female_mature,
sum(pop_female_total) sum_pop_female_total
from cb_popbystratum_sizegroup
group by survey_year
order by survey_year;

drop table cb_bioall_sizegroup;

create table cb_bioall_sizegroup as
select survey_year,
sum(bio_female_barren) sum_bio_female_barren,
sum(bio_female_trace) sum_bio_female_trace,
sum(bio_female_quarter) sum_bio_female_quarter,
sum(bio_female_half) sum_bio_female_half,
sum(bio_female_three_quarter) sum_bio_female_three_quarter,
sum(bio_female_full) sum_bio_female_full,
sum(bio_female_immature) sum_bio_female_immature,
sum(bio_female_mature) sum_bio_female_mature,
sum(bio_female_total) sum_bio_female_total
from cb_biobystratum_sizegroup
group by survey_year
order by survey_year;


drop table cb_varpop_sizegroup_sum;

create table cb_varpop_sizegroup_sum as
select distinct survey_year,
sum(varpop_female_barren) sum_varpop_female_barren,
sum(varpop_female_trace) sum_varpop_female_trace,
sum(varpop_female_quarter) sum_varpop_female_quarter,
sum(varpop_female_half) sum_varpop_female_half,
sum(varpop_female_three_quarter) sum_varpop_female_three_quarter,
sum(varpop_female_full) sum_varpop_female_full,
sum(varpop_female_immature) sum_varpop_female_immature,
sum(varpop_female_mature) sum_varpop_female_mature,
sum(varpop_female_total) sum_varpop_female_total
from cb_variancepop_sizegroup
group by survey_year
order by survey_year;

drop table cb_varbio_sizegroup_sum;

create table cb_varbio_sizegroup_sum as
select distinct survey_year,
sum(varbio_female_barren) sum_varbio_female_barren,
sum(varbio_female_trace) sum_varbio_female_trace,
sum(varbio_female_quarter) sum_varbio_female_quarter,
sum(varbio_female_half) sum_varbio_female_half,
sum(varbio_female_three_quarter) sum_varbio_female_three_quarter,
sum(varbio_female_full) sum_varbio_female_full,
sum(varbio_female_immature) sum_varbio_female_immature,
sum(varbio_female_mature) sum_varbio_female_mature,
sum(varbio_female_total) sum_varbio_female_total
from cb_variancebio_sizegroup
group by survey_year
order by survey_year;

drop table cb_pop_sizegroup_cv;

create table cb_pop_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_varpop_female_barren <> 0
	 then ((sqrt(sum_varpop_female_barren))/sum_pop_female_barren)
	 else 0
	 end) cv_pop_female_barren,
(CASE
	 when sum_pop_female_trace <> 0
	 then ((sqrt(sum_varpop_female_trace))/sum_pop_female_trace)
	 else 0
	 end) cv_pop_female_trace,   
(CASE
	 when sum_pop_female_quarter <> 0
	 then ((sqrt(sum_varpop_female_quarter))/sum_pop_female_quarter)
	 else 0
	 end) cv_pop_female_quarter,
(CASE
	 when sum_pop_female_half <> 0
	 then ((sqrt(sum_varpop_female_half))/sum_pop_female_half)
	 else 0
	 end) cv_pop_female_half,
(CASE
	 when sum_pop_female_three_quarter <> 0
	 then ((sqrt(sum_varpop_female_three_quarter))/sum_pop_female_three_quarter)
	 else 0
	 end) cv_pop_female_three_quarter,
(CASE
	 when sum_pop_female_full <> 0
	 then ((sqrt(sum_varpop_female_full))/sum_pop_female_full)
	 else 0
	 end) cv_pop_female_full,   
(CASE
	 when sum_pop_female_immature <> 0
	 then ((sqrt(sum_varpop_female_immature))/sum_pop_female_immature)
	 else 0
	 end) cv_pop_female_immature,	     
(CASE
	 when sum_pop_female_mature <> 0
	 then ((sqrt(sum_varpop_female_mature))/sum_pop_female_mature)
	 else 0
	 end) cv_pop_female_mature,	  
(CASE
	 when sum_pop_female_total <> 0
	 then ((sqrt(sum_varpop_female_total))/sum_pop_female_total)
	 else 0
	 end) cv_pop_female_total
     
from cb_varpop_sizegroup_sum a, cb_popall_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;	 	 	 

drop table cb_bio_sizegroup_cv;

create table cb_bio_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_bio_female_barren <> 0
	 then ((sqrt(sum_varbio_female_barren))/sum_bio_female_barren)
	 else 0
	 end) cv_bio_female_barren,
(CASE
	 when sum_bio_female_trace <> 0
	 then ((sqrt(sum_varbio_female_trace))/sum_bio_female_trace)
	 else 0
	 end) cv_bio_female_trace,   
(CASE
	 when sum_bio_female_quarter <> 0
	 then ((sqrt(sum_varbio_female_quarter))/sum_bio_female_quarter)
	 else 0
	 end) cv_bio_female_quarter,
(CASE
	 when sum_bio_female_half <> 0
	 then ((sqrt(sum_varbio_female_half))/sum_bio_female_half)
	 else 0
	 end) cv_bio_female_half,
(CASE
	 when sum_bio_female_three_quarter<> 0
	 then ((sqrt(sum_varbio_female_three_quarter))/sum_bio_female_three_quarter)
	 else 0
	 end) cv_bio_female_three_quarter,
(CASE
	 when sum_bio_female_full <> 0
	 then ((sqrt(sum_varbio_female_full))/sum_bio_female_full)
	 else 0
	 end) cv_bio_female_full,
(CASE
	 when sum_bio_female_immature <> 0
	 then ((sqrt(sum_varbio_female_immature))/sum_bio_female_immature)
	 else 0
	 end) cv_bio_female_immature,	     
(CASE
	 when sum_bio_female_mature <> 0
	 then ((sqrt(sum_varbio_female_mature))/sum_bio_female_mature)
	 else 0
	 end) cv_bio_female_mature,	  
(CASE
	 when sum_bio_female_total <> 0
	 then ((sqrt(sum_varbio_female_total))/sum_bio_female_total)
	 else 0
	 end) cv_bio_female_total  	 	 	 	  
from cb_varbio_sizegroup_sum a, cb_bioall_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;


-- CI calcs

create or replace view cb_sizegroup_stderr_pop as
select distinct survey_year,
(sqrt(sum_varpop_female_barren)) stderr_pop_female_barren,
(sqrt(sum_varpop_female_trace)) stderr_pop_female_trace,
(sqrt(sum_varpop_female_quarter)) stderr_pop_female_quarter,
(sqrt(sum_varpop_female_half)) stderr_pop_female_half,
(sqrt(sum_varpop_female_three_quarter)) stderr_pop_female_three_quarter,
(sqrt(sum_varpop_female_full)) stderr_pop_female_full,
(sqrt(sum_varpop_female_immature)) stderr_pop_female_immature,
(sqrt(sum_varpop_female_mature)) stderr_pop_female_mature,
(sqrt(sum_varpop_female_total)) stderr_pop_female_total
from cb_varpop_sizegroup_sum;

create or replace view cb_sizegroup_stderr_bio as
select distinct survey_year,
(sqrt(sum_varbio_female_barren)) stderr_bio_female_barren,
(sqrt(sum_varbio_female_trace)) stderr_bio_female_trace,
(sqrt(sum_varbio_female_quarter)) stderr_bio_female_quarter,
(sqrt(sum_varbio_female_half)) stderr_bio_female_half,
(sqrt(sum_varbio_female_three_quarter)) stderr_bio_female_three_quarter,
(sqrt(sum_varbio_female_full)) stderr_bio_female_full,
(sqrt(sum_varbio_female_immature)) stderr_bio_female_immature,
(sqrt(sum_varbio_female_mature)) stderr_bio_female_mature,
(sqrt(sum_varbio_female_total)) stderr_bio_female_total
from cb_varbio_sizegroup_sum;

drop table cb_sizegroup_confidence_pop;

create table cb_sizegroup_confidence_pop as
select distinct survey_year,
(1.96 * stderr_pop_female_barren) ci_pop_female_barren,
(1.96 * stderr_pop_female_trace) ci_pop_female_trace,
(1.96 * stderr_pop_female_quarter) ci_pop_female_quarter,
(1.96 * stderr_pop_female_half) ci_pop_female_half,
(1.96 * stderr_pop_female_three_quarter) ci_pop_female_three_quarter,
(1.96 * stderr_pop_female_full) ci_pop_female_full,
(1.96 * stderr_pop_female_immature) ci_pop_female_immature,
(1.96 * stderr_pop_female_mature) ci_pop_female_mature,
(1.96 * stderr_pop_female_total) ci_pop_female_total
from cb_sizegroup_stderr_pop;

drop table cb_sizegroup_confidence_bio;

create table cb_sizegroup_confidence_bio as
select distinct survey_year,
(1.96 * stderr_bio_female_barren) ci_bio_female_barren,
(1.96 * stderr_bio_female_trace) ci_bio_female_trace,
(1.96 * stderr_bio_female_quarter) ci_bio_female_quarter,
(1.96 * stderr_bio_female_half) ci_bio_female_half,
(1.96 * stderr_bio_female_three_quarter) ci_bio_female_three_quarter,
(1.96 * stderr_bio_female_full) ci_bio_female_full,
((1.96 * stderr_bio_female_immature)) ci_bio_female_immature,
((1.96 * stderr_bio_female_mature)) ci_bio_female_mature,
((1.96 * stderr_bio_female_total)) ci_bio_female_total
from cb_sizegroup_stderr_bio;



-- Final output for stocks

-- All Districts Combined

drop table cb_female_clutch_size_bio;

create table cb_female_clutch_size_bio as
select a.survey_year,
sum_bio_female_barren total_female_barren_biomass,cv_bio_female_barren cv_female_barren_mt, ci_bio_female_barren female_barren_biomass_ci,
sum_bio_female_trace total_female_trace_biomass,cv_bio_female_trace cv_female_trace_mt, ci_bio_female_trace female_trace_biomass_ci,
sum_bio_female_quarter total_female_quarter_biomass,cv_bio_female_quarter cv_female_quarter_mt, ci_bio_female_quarter female_quarter_biomass_ci,
sum_bio_female_half total_female_half_biomass,cv_bio_female_half cv_female_half_mt, ci_bio_female_half female_half_biomass_ci,
sum_bio_female_three_quarter total_female_three_quarter_biomass,cv_bio_female_three_quarter cv_female_three_quarter_mt, ci_bio_female_three_quarter female_three_quarter_biomass_ci,
sum_bio_female_full total_female_full_biomass,cv_bio_female_full cv_female_full_mt, ci_bio_female_full female_full_biomass_ci,
sum_bio_female_immature total_female_immature_biomass,cv_bio_female_immature cv_female_immature_mt, ci_bio_female_immature female_immature_biomass_ci,
sum_bio_female_mature total_female_mature_biomass,cv_bio_female_mature cv_female_mature_mt, ci_bio_female_mature female_mature_biomass_ci
from cb_bioall_sizegroup a, cb_sizegroup_confidence_bio b,cb_bio_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;


drop table cb_female_clutch_size_pop;

create table cb_female_clutch_size_pop as
select a.survey_year,
sum_pop_female_barren total_female_barren_popmass,cv_pop_female_barren cv_female_barren_mt, ci_pop_female_barren female_barren_popmass_ci,
sum_pop_female_trace total_female_trace_popmass,cv_pop_female_trace cv_female_trace_mt, ci_pop_female_trace female_trace_popmass_ci,
sum_pop_female_quarter total_female_quarter_popmass,cv_pop_female_quarter cv_female_quarter_mt, ci_pop_female_quarter female_quarter_popmass_ci,
sum_pop_female_half total_female_half_popmass,cv_pop_female_half cv_female_half_mt, ci_pop_female_half female_half_popmass_ci,
sum_pop_female_three_quarter total_female_three_quarter_popmass,cv_pop_female_three_quarter cv_female_three_quarter_mt, ci_pop_female_three_quarter female_three_quarter_popmass_ci,
sum_pop_female_full total_female_full_popmass,cv_pop_female_full cv_female_full_mt, ci_pop_female_full female_full_popmass_ci,
sum_pop_female_immature total_female_immature_popmass,cv_pop_female_immature cv_female_immature_mt, ci_pop_female_immature female_immature_popmass_ci,
sum_pop_female_mature total_female_mature_popmass,cv_pop_female_mature cv_female_mature_mt, ci_pop_female_mature female_mature_popmass_ci
from cb_popall_sizegroup a, cb_sizegroup_confidence_pop b,cb_pop_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;


alter table cb_female_clutch_size_pop
add SPECIES_CODE NUMBER;

update cb_female_clutch_size_pop
set SPECIES_CODE = 68560;

alter table cb_female_clutch_size_pop
add SURVEY_REGION VARCHAR2 (50 BYTE);

update cb_female_clutch_size_pop
set SURVEY_REGION = 'EBS';


alter table cb_female_clutch_size_bio
add SPECIES_CODE NUMBER;

update cb_female_clutch_size_bio
set SPECIES_CODE = 68560;

alter table cb_female_clutch_size_bio
add SURVEY_REGION VARCHAR2 (50 BYTE);

update cb_female_clutch_size_bio
set SURVEY_REGION = 'EBS';

