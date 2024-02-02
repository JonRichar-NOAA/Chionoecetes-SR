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
--create data table for immature bairdi
drop table immat_cbairdi;
create table immat_cbairdi as 
select * 
from crab.ebscrab
where sex = 1 and width <= 103
or sex = 2 and clutch_size = 0;


-- Don't want to average non-rkc in BB retow years, get rid of haul type 17 tows
drop table haul_newtimeseries_noretow;

create table haul_newtimeseries_noretow as
select * from haul_newtimeseries
where haul_type <> 17;


-- Create tables of raw catch by 1-mm size bin and sex
-- Separate by sex because male size group categories require shell condition
-- and female weights (post-2009) require clutch size

drop table cb_number_size1_male;

create table cb_number_size1_male as
select c.hauljoin,c.vessel,c.cruise,c.haul,h.gis_station,species_code,
shell_condition,(trunc(width/1) * 1)size1,
(sum(CASE
		 when species_code = 68560
		 and sex = 1
		 then sampling_factor
		 else 0
		 end)) number_male_size1
from crab.immat_cbairdi c, haul_newtimeseries_noretow h
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
		 (trunc(width/1) * 1);


-- Females (done separately from males because need clutch size info)

drop table cb_number_size1_female;

create table cb_number_size1_female as
select c.hauljoin,c.vessel,c.cruise,c.haul,h.gis_station,species_code,sex,clutch_size,
(trunc(width/1) * 1)size1,
(sum(CASE
		 when species_code = 68560
		 and sex = 2
		 then sampling_factor
		 else 0
		 end)) number_female_size1
from crab.immat_cbairdi c, haul_newtimeseries_noretow h
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
         sex,
         clutch_size,
		 (trunc(width/1) * 1);


-- unsexed

drop table cb_number_size1_unsexed;

create table cb_number_size1_unsexed as
select c.hauljoin,c.vessel,c.cruise,c.haul,h.gis_station,species_code,
(trunc(width/1) * 1)size1,
(sum(CASE
		 when species_code = 68560
		 and sex = 3
		 then sampling_factor
		 else 0
		 end)) number_unsexed_size1
from crab.immat_cbairdi c, haul_newtimeseries_noretow h
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
		 (trunc(width/1) * 1);



--  This section calculates the weight of the bairdi Tanner crab by haul, sex,
--  shell condition and 1-mm size group.  A width-weight regression
--  factor is applied, and multiplied by the number of crab caught in that
--  haul/sex/shellcon/size bin (from above section).  
--  The regression factor does not include unsexed crab, therefore no weights
--  will be calculated for unsexed crab


drop table cb_weight_grams_male;

create table cb_weight_grams_male as
select hauljoin,vessel,cruise,haul,gis_station,species_code,shell_condition,size1,
(CASE
--    WHEN cruise < 201001
--      THEN ((0.00019 * (power(size1,3.09894))) * number_male_size1)
    WHEN cruise >= 197501
      THEN ((0.00027 * (power(size1,3.022134))) * number_male_size1)
    ELSE 0
    END) wgt_male_size1
from cb_number_size1_male 
order by cruise,vessel,haul,gis_station,size1;

drop table cb_weight_grams_female;

create table cb_weight_grams_female as
select hauljoin,vessel,cruise,haul,gis_station,species_code,clutch_size,size1,
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
select hauljoin, vessel, cruise, haul, gis_station, species_code, size1,
	   (sum(CASE
	   			WHEN   sex = 2 
				   THEN number_female_size1
				ELSE 0
				END))  number_female_size1
                
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
				size1;
        
-- And calculate weight of crab by actual maturity        

drop table cb_weight_grams_matfem;

create table cb_weight_grams_matfem as
select hauljoin, vessel, cruise, haul, gis_station, species_code, size1,
	   (sum(CASE
	   			WHEN   clutch_size = 0  
				   THEN wgt_female_size1
				ELSE 0
				END))  wgt_female_size1	
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
  SIZE1                           NUMBER,
  SHELL_CONDITION                 NUMBER,
  WGT_MALE_SIZE1                  NUMBER,
  WGT_FEMALE_SIZE1                NUMBER
);
insert into cb_weight_grams_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,shell_condition,
wgt_male_size1,null
from cb_weight_grams_male;

insert into cb_weight_grams_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,
null,null,wgt_female_size1
from cb_weight_grams_matfem;

-- convert to metric tons

drop table cb_weight_mt_size1;

create table cb_weight_mt_size1 as
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,shell_condition,
(wgt_male_size1 * 0.000001) mt_male_size1,
(wgt_female_size1 * 0.000001) mt_female_size1
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
  SIZE1                           NUMBER,
  SHELL_CONDITION                 NUMBER,
  NUMBER_MALE_SIZE1               NUMBER,
  NUMBER_FEMALE_SIZE1            NUMBER,
  NUMBER_UNSEXED_SIZE1            NUMBER
);
insert into cb_number_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,shell_condition,
number_male_size1,null,null
from cb_number_size1_male;

insert into cb_number_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,null,
null,number_female_size1,null
from cb_number_size1_matfem;

insert into cb_number_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,null,
null,null,number_unsexed_size1
from cb_number_size1_unsexed;

-- This section sums the bairdi Tanner crab catch records by haul, sex,
-- and 1-mm size group.  

drop table cb_number_sizegroup;

create table cb_number_sizegroup as
select hauljoin, vessel, cruise, haul, gis_station, species_code, 
	   
  (sum(CASE
	   			WHEN  size1 between 0 and 29.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_le30,
        
  (sum(CASE
	   			WHEN  size1 between 0 and 39.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_le40,
    
    (sum(CASE
	   			WHEN  size1 between 0 and 49.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_le50,
        
    (sum(CASE
	   			WHEN  size1 between 0 and 59.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_le60,
        
    (sum(CASE
	   			WHEN  size1 between 29.9 and 39.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_30to40,
    (sum(CASE
	   			WHEN  size1 between 29.9 and 34.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_30to35,
    (sum(CASE
	   			WHEN  size1 between 34.9 and 39.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_35to40,          
    (sum(CASE
	   			WHEN  size1 between 39.9 and 49.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_40to50,
   (sum(CASE
	   			WHEN  size1 between 39.9 and 44.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_40to45,
    (sum(CASE
	   			WHEN  size1 between 44.9 and 49.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_45to50,
    (sum(CASE
	   			WHEN  size1 between 49.9 and 59.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_50to60,
        
    (sum(CASE
	   			WHEN  size1 between 29.9 and 49.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_30to50,
        
    (sum(CASE
	   			WHEN  size1 between 29.9 and 59.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_30to60,    
        
     (sum(CASE
	   			WHEN  size1 between 0 and 102.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_le103,
    (sum(CASE
	   			WHEN  size1 between 0 and 29.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_le30,     
    (sum(CASE
	   			WHEN  size1 between 0 and 39.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_le40,
    (sum(CASE
	   			WHEN  size1 between 0 and 49.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_le50,   
    (sum(CASE
	   			WHEN  size1 between 0 and 59.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_le60,   
    (sum(CASE
	   			WHEN  size1 between 29.9 and 39.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_30to40,
    (sum(CASE
	   			WHEN  size1 between 29.9 and 34.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_30to35,
    (sum(CASE
	   			WHEN  size1 between 34.9 and 39.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_35to40,
    (sum(CASE
	   			WHEN  size1 between 39.9 and 49.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_40to50,
    (sum(CASE
	   			WHEN  size1 between 39.9 and 44.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_40to45,
    (sum(CASE
	   			WHEN  size1 between 44.9 and 49.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_45to50,   
    (sum(CASE
	   			WHEN  size1 between 49.9 and 59.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_50to60,   
    (sum(CASE
	   			WHEN  size1 between 29.9 and 49.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_30to50,   
    (sum(CASE
	   			WHEN  size1 between 29.9 and 59.9
			    THEN number_female_size1
				ELSE 0
				END))  number_female_30to60, 
   sum(number_female_size1) number_female_immature
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
	   			WHEN  size1 between 0 and 29.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_le30,
    (sum(CASE
	   			WHEN  size1 between 0 and 39.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_le40,
    (sum(CASE
	   			WHEN  size1 between 0 and 49.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_le50,
    (sum(CASE
	   			WHEN  size1 between 0 and 59.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_le60,
    (sum(CASE
	   			WHEN  size1 between 29.9 and 39.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_30to40,
    (sum(CASE
	   			WHEN  size1 between 29.9 and 34.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_30to35,
    (sum(CASE
	   			WHEN  size1 between 34.9 and 39.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_35to40,
    (sum(CASE
	   			WHEN  size1 between 39.9 and 49.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_40to50,
    (sum(CASE
	   			WHEN  size1 between 39.9 and 44.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_40to45,
    (sum(CASE
	   			WHEN  size1 between 44.9 and 49.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_45to50,
   (sum(CASE
	   			WHEN  size1 between 49.9 and 59.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_50to60,
    (sum(CASE
	   			WHEN  size1 between 29.9 and 49.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_30to50,  
     (sum(CASE
	   			WHEN  size1 between 29.9 and 59.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_30to60, 
    (sum(CASE
	   			WHEN  size1 between 0.0 and 102.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_le103,              
    (sum(CASE
	   			WHEN  size1 between 0 and 29.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_le30,
    (sum(CASE
	   			WHEN  size1 between 0 and 39.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_le40,
    (sum(CASE
	   			WHEN  size1 between 0 and 49.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_le50,
    (sum(CASE
	   			WHEN  size1 between 0 and 59.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_le60,
    (sum(CASE
	   			WHEN  size1 between 29.9 and 39.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_30to40,
    (sum(CASE
	   			WHEN  size1 between 29.9 and 34.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_30to35,
    (sum(CASE
	   			WHEN  size1 between 34.9 and 39.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_35to40,
    (sum(CASE
	   			WHEN  size1 between 39.9 and 49.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_40to50,
    (sum(CASE
	   			WHEN  size1 between 39.9 and 44.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_40to45,
    (sum(CASE
	   			WHEN  size1 between 44.9 and 49.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_45to50,
   (sum(CASE
	   			WHEN  size1 between 49.9 and 59.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_50to60,
    (sum(CASE
	   			WHEN  size1 between 29.9 and 49.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_30to50,  
     (sum(CASE
	   			WHEN  size1 between 29.9 and 59.9
			    THEN mt_female_size1
				ELSE 0
				END))  mt_female_30to60,          
     sum(mt_female_size1) mt_female_immature
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
select h.hauljoin,h.vessel,h.cruise,h.haul,h.gis_station,survey_year,
nvl(species_code,68560) species_code,
nvl(number_male_le30,0) number_male_le30,
nvl(number_male_le40,0) number_male_le40,
nvl(number_male_le50,0) number_male_le50,
nvl(number_male_le60,0) number_male_le60,
nvl(number_male_30to40,0) number_male_30to40,
nvl(number_male_30to35,0) number_male_30to35,
nvl(number_male_35to40,0) number_male_35to40,
nvl(number_male_40to50,0) number_male_40to50,
nvl(number_male_40to45,0) number_male_40to45,
nvl(number_male_45to50,0) number_male_45to50,
nvl(number_male_50to60,0) number_male_50to60,
nvl(number_male_30to50,0) number_male_30to50,
nvl(number_male_30to60,0) number_male_30to60,
nvl(number_male_le103,0) number_male_le103,
nvl(number_female_le30,0) number_female_le30,
nvl(number_female_le40,0) number_female_le40,
nvl(number_female_le50,0) number_female_le50,
nvl(number_female_le60,0) number_female_le60,
nvl(number_female_30to40,0) number_female_30to40,
nvl(number_female_30to35,0) number_female_30to35,
nvl(number_female_35to40,0) number_female_35to40,
nvl(number_female_40to50,0) number_female_40to50,
nvl(number_female_40to45,0) number_female_40to45,
nvl(number_female_45to50,0) number_female_45to50,
nvl(number_female_50to60,0) number_female_50to60,
nvl(number_female_30to50,0) number_female_30to50,
nvl(number_female_30to60,0) number_female_30to60,
nvl(number_female_immature,0) number_female_immature
from haul_newtimeseries_noretow h full outer join cb_number_sizegroup c
on h.hauljoin = c.hauljoin
where haul_type <> 17;

--  Similarly, by weight.

drop table cb_wgt_sizegroup_union;

create table cb_wgt_sizegroup_union as
select h.hauljoin,h.vessel,h.cruise,h.haul,h.gis_station,survey_year,
nvl(species_code,68560) species_code,
nvl(mt_male_le30,0) mt_male_le30,
nvl(mt_male_le40,0) mt_male_le40,
nvl(mt_male_le50,0) mt_male_le50,
nvl(mt_male_le60,0) mt_male_le60,
nvl(mt_male_30to40,0) mt_male_30to40,
nvl(mt_male_30to35,0) mt_male_30to35,
nvl(mt_male_35to40,0) mt_male_35to40,
nvl(mt_male_40to50,0) mt_male_40to50,
nvl(mt_male_40to45,0) mt_male_40to45,
nvl(mt_male_45to50,0) mt_male_45to50,
nvl(mt_male_50to60,0) mt_male_50to60,
nvl(mt_male_30to50,0) mt_male_30to50,
nvl(mt_male_30to60,0) mt_male_30to60,
nvl(mt_male_le103,0) mt_male_le103,
nvl(mt_female_le30,0) mt_female_le30,
nvl(mt_female_le40,0) mt_female_le40,
nvl(mt_female_le50,0) mt_female_le50,
nvl(mt_female_le60,0) mt_female_le60,
nvl(mt_female_30to40,0) mt_female_30to40,
nvl(mt_female_30to35,0) mt_female_30to35,
nvl(mt_female_35to40,0) mt_female_35to40,
nvl(mt_female_40to50,0) mt_female_40to50,
nvl(mt_female_40to45,0) mt_female_40to45,
nvl(mt_female_45to50,0) mt_female_45to50,
nvl(mt_female_50to60,0) mt_female_50to60,
nvl(mt_female_30to50,0) mt_female_30to50,
nvl(mt_female_30to60,0) mt_female_30to60,
nvl(mt_female_immature,0) mt_female_immature
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
select c.hauljoin,c.vessel,c.cruise,c.haul,c.gis_station,c.survey_year,c.species_code,
(number_male_le30 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_le30,
(number_male_le40 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_le40,
(number_male_le50 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_le50,
(number_male_le60 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_le60,
(number_male_30to40 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_30to40,
(number_male_30to35 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_30to35,
(number_male_35to40 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_35to40,
(number_male_40to50 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_40to50,
(number_male_40to45 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_40to45,
(number_male_45to50 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_45to50,
(number_male_50to60 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_50to60,
(number_male_30to50 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_30to50,
(number_male_30to60 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_30to60,
(number_male_le103 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_le103,
(number_female_le30 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_le30,
(number_female_le40 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_le40,
(number_female_le50 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_le50,
(number_female_le60 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_le60,
(number_female_30to40 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_30to40,
(number_female_30to35 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_30to35,
(number_female_35to40 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_35to40,
(number_female_40to50 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_40to50,
(number_female_40to45 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_40to45,
(number_female_45to50 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_45to50,
(number_female_50to60 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_50to60,
(number_female_30to50 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_30to50,
(number_female_30to60 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_30to60,
(number_female_immature / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_immature
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
select c.hauljoin,c.vessel,c.cruise,c.haul,mid_latitude,mid_longitude,
c.gis_station,c.survey_year,c.species_code,
(mt_male_le30 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_le30,
(mt_male_le40 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_le40,
(mt_male_le50 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_le50,
(mt_male_le60 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_le60,
(mt_male_30to40 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_30to40,
(mt_male_30to35 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_30to35,
(mt_male_35to40 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_35to40,
(mt_male_40to50 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_40to50,
(mt_male_40to45 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_40to45,
(mt_male_45to50 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_45to50,
(mt_male_50to60 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_50to60,
(mt_male_30to50 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_30to50,
(mt_male_30to60 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_30to60,
(mt_male_le103 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_le103,
(mt_female_le30 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_le30,
(mt_female_le40 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_le40,
(mt_female_le50 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_le50,
(mt_female_le60 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_le60,
(mt_female_30to40 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_30to40,
(mt_female_30to35 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_30to35,
(mt_female_35to40 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_35to40,
(mt_female_40to50 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_40to50,
(mt_female_40to45 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_40to45,
(mt_female_45to50 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_45to50,
(mt_female_50to60 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_50to60,
(mt_female_30to50 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_30to50,
(mt_female_30to60 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_30to60,
(mt_female_immature / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_immature
from cb_wgt_sizegroup_union c, haul_newtimeseries_noretow h
where c.hauljoin = h.hauljoin
and haul_type <> 17;


drop table cb_meancpuenum_sizegroup;

create table cb_meancpuenum_sizegroup as
select c.survey_year,district,
AVG (male_cpuenum_le30) meancpuenum_male_le30,
AVG (male_cpuenum_le40) meancpuenum_male_le40,
AVG (male_cpuenum_le50) meancpuenum_male_le50,
AVG (male_cpuenum_le60) meancpuenum_male_le60,
AVG (male_cpuenum_30to40) meancpuenum_male_30to40,
AVG (male_cpuenum_30to35) meancpuenum_male_30to35,
AVG (male_cpuenum_35to40) meancpuenum_male_35to40,
AVG (male_cpuenum_40to50) meancpuenum_male_40to50,
AVG (male_cpuenum_40to45) meancpuenum_male_40to45,
AVG (male_cpuenum_45to50) meancpuenum_male_45to50,
AVG (male_cpuenum_50to60) meancpuenum_male_50to60,
AVG (male_cpuenum_30to50) meancpuenum_male_30to50,
AVG (male_cpuenum_30to60) meancpuenum_male_30to60,
AVG (male_cpuenum_le103) meancpuenum_male_le103,
AVG (female_cpuenum_le30) meancpuenum_female_le30,
AVG (female_cpuenum_le40) meancpuenum_female_le40,
AVG (female_cpuenum_le50) meancpuenum_female_le50,
AVG (female_cpuenum_le60) meancpuenum_female_le60,
AVG (female_cpuenum_30to40) meancpuenum_female_30to40,
AVG (female_cpuenum_30to35) meancpuenum_female_30to35,
AVG (female_cpuenum_35to40) meancpuenum_female_35to40,
AVG (female_cpuenum_40to50) meancpuenum_female_40to50,
AVG (female_cpuenum_40to45) meancpuenum_female_40to45,
AVG (female_cpuenum_45to50) meancpuenum_female_45to50,
AVG (female_cpuenum_50to60) meancpuenum_female_50to60,
AVG (female_cpuenum_30to50) meancpuenum_female_30to50,
AVG (female_cpuenum_30to60) meancpuenum_female_30to60,
AVG (female_cpuenum_immature) meancpuenum_female_immature
from cb_cpuenum_sizegroup c, strata_bairdi_newtimeseries s
where c.gis_station = s.station_id
and c.survey_year = s.survey_year
group by c.survey_year,district;

drop table cb_meancpuewgt_sizegroup;

create table cb_meancpuewgt_sizegroup as
select c.survey_year,district,
AVG (male_cpuewgt_le30) meancpuewgt_male_le30,
AVG (male_cpuewgt_le40) meancpuewgt_male_le40,
AVG (male_cpuewgt_le50) meancpuewgt_male_le50,
AVG (male_cpuewgt_le60) meancpuewgt_male_le60,
AVG (male_cpuewgt_30to40) meancpuewgt_male_30to40,
AVG (male_cpuewgt_30to35) meancpuewgt_male_30to35,
AVG (male_cpuewgt_35to40) meancpuewgt_male_35to40,
AVG (male_cpuewgt_40to50) meancpuewgt_male_40to50,
AVG (male_cpuewgt_40to45) meancpuewgt_male_40to45,
AVG (male_cpuewgt_45to50) meancpuewgt_male_45to50,
AVG (male_cpuewgt_50to60) meancpuewgt_male_50to60,
AVG (male_cpuewgt_30to50) meancpuewgt_male_30to50,
AVG (male_cpuewgt_30to60) meancpuewgt_male_30to60,
AVG (male_cpuewgt_le103) meancpuewgt_male_le103,
AVG (female_cpuewgt_le30) meancpuewgt_female_le30,
AVG (female_cpuewgt_le40) meancpuewgt_female_le40,
AVG (female_cpuewgt_le50) meancpuewgt_female_le50,
AVG (female_cpuewgt_le60) meancpuewgt_female_le60,
AVG (female_cpuewgt_30to40) meancpuewgt_female_30to40,
AVG (female_cpuewgt_30to35) meancpuewgt_female_30to35,
AVG (female_cpuewgt_35to40) meancpuewgt_female_35to40,
AVG (female_cpuewgt_40to50) meancpuewgt_female_40to50,
AVG (female_cpuewgt_40to45) meancpuewgt_female_40to45,
AVG (female_cpuewgt_45to50) meancpuewgt_female_45to50,
AVG (female_cpuewgt_50to60) meancpuewgt_female_50to60,
AVG (female_cpuewgt_30to50) meancpuewgt_female_30to50,
AVG (female_cpuewgt_30to60) meancpuewgt_female_30to60,
AVG (female_cpuewgt_immature) meancpuewgt_female_immature
from cb_cpuewgt_sizegroup c, strata_bairdi_newtimeseries s
where c.gis_station = s.station_id
and c.survey_year = s.survey_year
group by c.survey_year,district;


drop table cb_popbystratum_sizegroup;

create table cb_popbystratum_sizegroup as
select distinct c.survey_year,stratum,c.district,
 (meancpuenum_male_le30 * total_area) pop_male_le30,
 (meancpuenum_male_le40 * total_area) pop_male_le40,
 (meancpuenum_male_le50 * total_area) pop_male_le50,
 (meancpuenum_male_le60 * total_area) pop_male_le60,
 (meancpuenum_male_30to40 * total_area) pop_male_30to40,
 (meancpuenum_male_30to35 * total_area) pop_male_30to35,
 (meancpuenum_male_35to40 * total_area) pop_male_35to40,
 (meancpuenum_male_40to50 * total_area) pop_male_40to50,
 (meancpuenum_male_40to45 * total_area) pop_male_40to45,
 (meancpuenum_male_45to50 * total_area) pop_male_45to50,
 (meancpuenum_male_50to60 * total_area) pop_male_50to60,
 (meancpuenum_male_30to50 * total_area) pop_male_30to50,
 (meancpuenum_male_30to60 * total_area) pop_male_30to60,
 (meancpuenum_male_le103 * total_area) pop_male_le103,
 (meancpuenum_female_le30 * total_area) pop_female_le30,
 (meancpuenum_female_le40 * total_area) pop_female_le40,
 (meancpuenum_female_le50 * total_area) pop_female_le50,
 (meancpuenum_female_le60 * total_area) pop_female_le60,
 (meancpuenum_female_30to40 * total_area) pop_female_30to40,
 (meancpuenum_female_30to35 * total_area) pop_female_30to35,
 (meancpuenum_female_35to40 * total_area) pop_female_35to40,
 (meancpuenum_female_40to50 * total_area) pop_female_40to50,
 (meancpuenum_female_40to45 * total_area) pop_female_40to45,
 (meancpuenum_female_45to50 * total_area) pop_female_45to50,
 (meancpuenum_female_50to60 * total_area) pop_female_50to60,
 (meancpuenum_female_30to50 * total_area) pop_female_30to50,
 (meancpuenum_female_30to60 * total_area) pop_female_30to60,
 (meancpuenum_female_immature * total_area) pop_female_immature
from cb_meancpuenum_sizegroup c, strata_bairdi_newtimeseries s
where c.district = s.district
and c.survey_year = s.survey_year
order by survey_year,district;

drop table cb_biobystratum_sizegroup;

create table cb_biobystratum_sizegroup as
select distinct c.survey_year,stratum,c.district,
(meancpuewgt_male_le30 * total_area) bio_male_le30,
(meancpuewgt_male_le40 * total_area) bio_male_le40,
(meancpuewgt_male_le50 * total_area) bio_male_le50,
(meancpuewgt_male_le60 * total_area) bio_male_le60,
(meancpuewgt_male_30to40 * total_area) bio_male_30to40,
(meancpuewgt_male_30to35 * total_area) bio_male_30to35,
(meancpuewgt_male_35to40 * total_area) bio_male_35to40,
(meancpuewgt_male_40to50 * total_area) bio_male_40to50,
(meancpuewgt_male_40to45 * total_area) bio_male_40to45,
(meancpuewgt_male_45to50 * total_area) bio_male_45to50,
(meancpuewgt_male_50to60 * total_area) bio_male_50to60,
(meancpuewgt_male_30to50 * total_area) bio_male_30to50,
(meancpuewgt_male_30to60 * total_area) bio_male_30to60,
(meancpuewgt_male_le103 * total_area) bio_male_le103,
(meancpuewgt_female_le30 * total_area) bio_female_le30,
(meancpuewgt_female_le40 * total_area) bio_female_le40,
(meancpuewgt_female_le50 * total_area) bio_female_le50,
(meancpuewgt_female_le60 * total_area) bio_female_le60,
(meancpuewgt_female_30to40 * total_area) bio_female_30to40,
(meancpuewgt_female_30to35 * total_area) bio_female_30to35,
(meancpuewgt_female_35to40 * total_area) bio_female_35to40,
(meancpuewgt_female_40to50 * total_area) bio_female_40to50,
(meancpuewgt_female_40to45 * total_area) bio_female_40to45,
(meancpuewgt_female_45to50 * total_area) bio_female_45to50,
(meancpuewgt_female_50to60 * total_area) bio_female_50to60,
(meancpuewgt_female_30to50 * total_area) bio_female_30to50,
(meancpuewgt_female_30to60 * total_area) bio_female_30to60,
(meancpuewgt_female_immature * total_area) bio_female_immature
from cb_meancpuewgt_sizegroup c, strata_bairdi_newtimeseries s
where c.district = s.district
and c.survey_year = s.survey_year
order by survey_year,district;

drop table cb_varcpuenum_sizegroup;

create table cb_varcpuenum_sizegroup as
select c.survey_year,district,
VARIANCE (male_cpuenum_le30) varcpuenum_male_le30,
VARIANCE (male_cpuenum_le40) varcpuenum_male_le40,
VARIANCE (male_cpuenum_le50) varcpuenum_male_le50,
VARIANCE (male_cpuenum_le60) varcpuenum_male_le60,
VARIANCE (male_cpuenum_30to40) varcpuenum_male_30to40,
VARIANCE (male_cpuenum_30to35) varcpuenum_male_30to35,
VARIANCE (male_cpuenum_35to40) varcpuenum_male_35to40,
VARIANCE (male_cpuenum_40to50) varcpuenum_male_40to50,
VARIANCE (male_cpuenum_40to45) varcpuenum_male_40to45,
VARIANCE (male_cpuenum_45to50) varcpuenum_male_45to50,
VARIANCE (male_cpuenum_50to60) varcpuenum_male_50to60,
VARIANCE (male_cpuenum_30to50) varcpuenum_male_30to50,
VARIANCE (male_cpuenum_30to60) varcpuenum_male_30to60,
VARIANCE (male_cpuenum_le103) varcpuenum_male_le103,
VARIANCE (female_cpuenum_le30) varcpuenum_female_le30,
VARIANCE (female_cpuenum_le40) varcpuenum_female_le40,
VARIANCE (female_cpuenum_le50) varcpuenum_female_le50,
VARIANCE (female_cpuenum_le60) varcpuenum_female_le60,
VARIANCE (female_cpuenum_30to40) varcpuenum_female_30to40,
VARIANCE (female_cpuenum_30to35) varcpuenum_female_30to35,
VARIANCE (female_cpuenum_35to40) varcpuenum_female_35to40,
VARIANCE (female_cpuenum_40to50) varcpuenum_female_40to50,
VARIANCE (female_cpuenum_40to45) varcpuenum_female_40to45,
VARIANCE (female_cpuenum_45to50) varcpuenum_female_45to50,
VARIANCE (female_cpuenum_50to60) varcpuenum_female_50to60,
VARIANCE (female_cpuenum_30to50) varcpuenum_female_30to50,
VARIANCE (female_cpuenum_30to60) varcpuenum_female_30to60,
VARIANCE (female_cpuenum_immature) varcpuenum_female_immature
from cb_cpuenum_sizegroup c, strata_bairdi_newtimeseries s
where c.gis_station = s.station_id
and c.survey_year = s.survey_year
group by c.survey_year,district;

drop table cb_varcpuewgt_sizegroup;

create table cb_varcpuewgt_sizegroup as
select c.survey_year,district,
VARIANCE (male_cpuewgt_le30) varcpuewgt_male_le30,
VARIANCE (male_cpuewgt_le40) varcpuewgt_male_le40,
VARIANCE (male_cpuewgt_le50) varcpuewgt_male_le50,
VARIANCE (male_cpuewgt_le60) varcpuewgt_male_le60,
VARIANCE (male_cpuewgt_30to40) varcpuewgt_male_30to40,
VARIANCE (male_cpuewgt_30to35) varcpuewgt_male_30to35,
VARIANCE (male_cpuewgt_35to40) varcpuewgt_male_35to40,
VARIANCE (male_cpuewgt_40to50) varcpuewgt_male_40to50,
VARIANCE (male_cpuewgt_40to45) varcpuewgt_male_40to45,
VARIANCE (male_cpuewgt_45to50) varcpuewgt_male_45to50,
VARIANCE (male_cpuewgt_50to60) varcpuewgt_male_50to60,
VARIANCE (male_cpuewgt_30to50) varcpuewgt_male_30to50,
VARIANCE (male_cpuewgt_30to60) varcpuewgt_male_30to60,
VARIANCE (male_cpuewgt_le103) varcpuewgt_male_le103,
VARIANCE (female_cpuewgt_le30) varcpuewgt_female_le30,
VARIANCE (female_cpuewgt_le40) varcpuewgt_female_le40,
VARIANCE (female_cpuewgt_le50) varcpuewgt_female_le50,
VARIANCE (female_cpuewgt_le60) varcpuewgt_female_le60,
VARIANCE (female_cpuewgt_30to40) varcpuewgt_female_30to40,
VARIANCE (female_cpuewgt_30to35) varcpuewgt_female_30to35,
VARIANCE (female_cpuewgt_35to40) varcpuewgt_female_35to40,
VARIANCE (female_cpuewgt_40to50) varcpuewgt_female_40to50,
VARIANCE (female_cpuewgt_40to45) varcpuewgt_female_40to45,
VARIANCE (female_cpuewgt_45to50) varcpuewgt_female_45to50,
VARIANCE (female_cpuewgt_50to60) varcpuewgt_female_50to60,
VARIANCE (female_cpuewgt_30to50) varcpuewgt_female_30to50,
VARIANCE (female_cpuewgt_30to60) varcpuewgt_female_30to60,
VARIANCE (female_cpuewgt_immature) varcpuewgt_female_immature
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
((varcpuenum_male_le30 * (power(total_area,2)))/number_tows) varpop_male_le30,
((varcpuenum_male_le40 * (power(total_area,2)))/number_tows) varpop_male_le40,
((varcpuenum_male_le50 * (power(total_area,2)))/number_tows) varpop_male_le50,
((varcpuenum_male_le60 * (power(total_area,2)))/number_tows) varpop_male_le60,
((varcpuenum_male_30to40 * (power(total_area,2)))/number_tows) varpop_male_30to40,
((varcpuenum_male_30to35 * (power(total_area,2)))/number_tows) varpop_male_30to35,
((varcpuenum_male_35to40 * (power(total_area,2)))/number_tows) varpop_male_35to40,
((varcpuenum_male_40to50 * (power(total_area,2)))/number_tows) varpop_male_40to50,
((varcpuenum_male_40to45 * (power(total_area,2)))/number_tows) varpop_male_40to45,
((varcpuenum_male_45to50 * (power(total_area,2)))/number_tows) varpop_male_45to50,
((varcpuenum_male_50to60 * (power(total_area,2)))/number_tows) varpop_male_50to60,
((varcpuenum_male_30to50 * (power(total_area,2)))/number_tows) varpop_male_30to50,
((varcpuenum_male_30to60 * (power(total_area,2)))/number_tows) varpop_male_30to60,
((varcpuenum_male_le103 * (power(total_area,2)))/number_tows) varpop_male_le103,
((varcpuenum_female_le30 * (power(total_area,2)))/number_tows) varpop_female_le30,
((varcpuenum_female_le40 * (power(total_area,2)))/number_tows) varpop_female_le40,
((varcpuenum_female_le50 * (power(total_area,2)))/number_tows) varpop_female_le50,
((varcpuenum_female_le60 * (power(total_area,2)))/number_tows) varpop_female_le60,
((varcpuenum_female_30to40 * (power(total_area,2)))/number_tows) varpop_female_30to40,
((varcpuenum_female_30to35 * (power(total_area,2)))/number_tows) varpop_female_30to35,
((varcpuenum_female_35to40 * (power(total_area,2)))/number_tows) varpop_female_35to40,
((varcpuenum_female_40to50 * (power(total_area,2)))/number_tows) varpop_female_40to50,
((varcpuenum_female_40to45 * (power(total_area,2)))/number_tows) varpop_female_40to45,
((varcpuenum_female_45to50 * (power(total_area,2)))/number_tows) varpop_female_45to50,
((varcpuenum_female_50to60 * (power(total_area,2)))/number_tows) varpop_female_50to60,
((varcpuenum_female_30to50 * (power(total_area,2)))/number_tows) varpop_female_30to50,
((varcpuenum_female_30to60 * (power(total_area,2)))/number_tows) varpop_female_30to60,
((varcpuenum_female_immature * (power(total_area,2)))/number_tows) varpop_female_immature
from strata_bairdi_newtimeseries s, cb_varcpuenum_sizegroup c, cb_haulcount n
where c.district = s.district
and c.district = n.district
and c.survey_year = s.survey_year
and c.survey_year = n.survey_year
order by c.survey_year,stratum;

drop table cb_variancebio_sizegroup;

create table cb_variancebio_sizegroup as
select distinct c.survey_year,stratum,c.district,
((varcpuewgt_male_le30 * (power(total_area,2)))/number_tows) varbio_male_le30,
((varcpuewgt_male_le40 * (power(total_area,2)))/number_tows) varbio_male_le40,
((varcpuewgt_male_le50 * (power(total_area,2)))/number_tows) varbio_male_le50,
((varcpuewgt_male_le60 * (power(total_area,2)))/number_tows) varbio_male_le60,
((varcpuewgt_male_30to40 * (power(total_area,2)))/number_tows) varbio_male_30to40,
((varcpuewgt_male_30to35 * (power(total_area,2)))/number_tows) varbio_male_30to35,
((varcpuewgt_male_35to40 * (power(total_area,2)))/number_tows) varbio_male_35to40,
((varcpuewgt_male_40to50 * (power(total_area,2)))/number_tows) varbio_male_40to50,
((varcpuewgt_male_40to45 * (power(total_area,2)))/number_tows) varbio_male_40to45,
((varcpuewgt_male_45to50 * (power(total_area,2)))/number_tows) varbio_male_45to50,
((varcpuewgt_male_50to60 * (power(total_area,2)))/number_tows) varbio_male_50to60,
((varcpuewgt_male_30to50 * (power(total_area,2)))/number_tows) varbio_male_30to50,
((varcpuewgt_male_30to60 * (power(total_area,2)))/number_tows) varbio_male_30to60,
((varcpuewgt_male_le103 * (power(total_area,2)))/number_tows) varbio_male_le103,
((varcpuewgt_female_le30 * (power(total_area,2)))/number_tows) varbio_female_le30,
((varcpuewgt_female_le40 * (power(total_area,2)))/number_tows) varbio_female_le40,
((varcpuewgt_female_le50 * (power(total_area,2)))/number_tows) varbio_female_le50,
((varcpuewgt_female_le60 * (power(total_area,2)))/number_tows) varbio_female_le60,
((varcpuewgt_female_30to40 * (power(total_area,2)))/number_tows) varbio_female_30to40,
((varcpuewgt_female_30to35 * (power(total_area,2)))/number_tows) varbio_female_30to35,
((varcpuewgt_female_35to40 * (power(total_area,2)))/number_tows) varbio_female_35to40,
((varcpuewgt_female_40to50 * (power(total_area,2)))/number_tows) varbio_female_40to50,
((varcpuewgt_female_40to45 * (power(total_area,2)))/number_tows) varbio_female_40to45,
((varcpuewgt_female_45to50 * (power(total_area,2)))/number_tows) varbio_female_45to50,
((varcpuewgt_female_50to60 * (power(total_area,2)))/number_tows) varbio_female_50to60,
((varcpuewgt_female_30to50 * (power(total_area,2)))/number_tows) varbio_female_30to50,
((varcpuewgt_female_30to60 * (power(total_area,2)))/number_tows) varbio_female_30to60,
((varcpuewgt_female_immature * (power(total_area,2)))/number_tows) varbio_female_immature
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
sum(pop_male_le30) sum_pop_male_le30,
sum(pop_male_le40) sum_pop_male_le40,
sum(pop_male_le50) sum_pop_male_le50,
sum(pop_male_le60) sum_pop_male_le60,
sum(pop_male_30to40) sum_pop_male_30to40,
sum(pop_male_30to35) sum_pop_male_30to35,
sum(pop_male_35to40) sum_pop_male_35to40,
sum(pop_male_40to50) sum_pop_male_40to50,
sum(pop_male_40to45) sum_pop_male_40to45,
sum(pop_male_45to50) sum_pop_male_45to50,
sum(pop_male_50to60) sum_pop_male_50to60,
sum(pop_male_30to50) sum_pop_male_30to50,
sum(pop_male_30to60) sum_pop_male_30to60,
sum(pop_male_le103) sum_pop_male_le103,
sum(pop_female_le30) sum_pop_female_le30,
sum(pop_female_le40) sum_pop_female_le40,
sum(pop_female_le50) sum_pop_female_le50,
sum(pop_female_le60) sum_pop_female_le60,
sum(pop_female_30to40) sum_pop_female_30to40,
sum(pop_female_30to35) sum_pop_female_30to35,
sum(pop_female_35to40) sum_pop_female_35to40,
sum(pop_female_40to50) sum_pop_female_40to50,
sum(pop_female_40to45) sum_pop_female_40to45,
sum(pop_female_45to50) sum_pop_female_45to50,
sum(pop_female_50to60) sum_pop_female_50to60,
sum(pop_female_30to50) sum_pop_female_30to50,
sum(pop_female_30to60) sum_pop_female_30to60,
sum(pop_female_immature) sum_pop_female_immature
from cb_popbystratum_sizegroup
group by survey_year
order by survey_year;

drop table cb_bioall_sizegroup;

create table cb_bioall_sizegroup as
select survey_year,
sum(bio_male_le30) sum_bio_male_le30,
sum(bio_male_le40) sum_bio_male_le40,
sum(bio_male_le50) sum_bio_male_le50,
sum(bio_male_le60) sum_bio_male_le60,
sum(bio_male_30to40) sum_bio_male_30to40,
sum(bio_male_30to35) sum_bio_male_30to35,
sum(bio_male_35to40) sum_bio_male_35to40,
sum(bio_male_40to50) sum_bio_male_40to50,
sum(bio_male_40to45) sum_bio_male_40to45,
sum(bio_male_45to50) sum_bio_male_45to50,
sum(bio_male_50to60) sum_bio_male_50to60,
sum(bio_male_30to50) sum_bio_male_30to50,
sum(bio_male_30to60) sum_bio_male_30to60,
sum(bio_male_le103) sum_bio_male_le103,
sum(bio_female_le30) sum_bio_female_le30,
sum(bio_female_le40) sum_bio_female_le40,
sum(bio_female_le50) sum_bio_female_le50,
sum(bio_female_le60) sum_bio_female_le60,
sum(bio_female_30to40) sum_bio_female_30to40,
sum(bio_female_30to35) sum_bio_female_30to35,
sum(bio_female_35to40) sum_bio_female_35to40,
sum(bio_female_40to50) sum_bio_female_40to50,
sum(bio_female_40to45) sum_bio_female_40to45,
sum(bio_female_45to50) sum_bio_female_45to50,
sum(bio_female_50to60) sum_bio_female_50to60,
sum(bio_female_30to50) sum_bio_female_30to50,
sum(bio_female_30to60) sum_bio_female_30to60,
sum(bio_female_immature) sum_bio_female_immature
from cb_biobystratum_sizegroup
group by survey_year
order by survey_year;


drop table cb_varpop_sizegroup_sum;

create table cb_varpop_sizegroup_sum as
select distinct survey_year,
sum(varpop_male_le30) sum_varpop_male_le30,
sum(varpop_male_le40) sum_varpop_male_le40,
sum(varpop_male_le50) sum_varpop_male_le50,
sum(varpop_male_le60) sum_varpop_male_le60,
sum(varpop_male_30to40) sum_varpop_male_30to40,
sum(varpop_male_30to35) sum_varpop_male_30to35,
sum(varpop_male_35to40) sum_varpop_male_35to40,
sum(varpop_male_40to50) sum_varpop_male_40to50,
sum(varpop_male_40to45) sum_varpop_male_40to45,
sum(varpop_male_45to50) sum_varpop_male_45to50,
sum(varpop_male_50to60) sum_varpop_male_50to60,
sum(varpop_male_30to50) sum_varpop_male_30to50,
sum(varpop_male_30to60) sum_varpop_male_30to60,
sum(varpop_male_le103) sum_varpop_male_le103,
sum(varpop_female_le30) sum_varpop_female_le30,
sum(varpop_female_le40) sum_varpop_female_le40,
sum(varpop_female_le50) sum_varpop_female_le50,
sum(varpop_female_le60) sum_varpop_female_le60,
sum(varpop_female_30to40) sum_varpop_female_30to40,
sum(varpop_female_30to35) sum_varpop_female_30to35,
sum(varpop_female_35to40) sum_varpop_female_35to40,
sum(varpop_female_40to50) sum_varpop_female_40to50,
sum(varpop_female_40to45) sum_varpop_female_40to45,
sum(varpop_female_45to50) sum_varpop_female_45to50,
sum(varpop_female_50to60) sum_varpop_female_50to60,
sum(varpop_female_30to50) sum_varpop_female_30to50,
sum(varpop_female_30to60) sum_varpop_female_30to60,
sum(varpop_female_immature) sum_varpop_female_immature
from cb_variancepop_sizegroup
group by survey_year
order by survey_year;

drop table cb_varbio_sizegroup_sum;

create table cb_varbio_sizegroup_sum as
select distinct survey_year,
sum(varbio_male_le30) sum_varbio_male_le30,
sum(varbio_male_le40) sum_varbio_male_le40,
sum(varbio_male_le50) sum_varbio_male_le50,
sum(varbio_male_le60) sum_varbio_male_le60,
sum(varbio_male_30to40) sum_varbio_male_30to40,
sum(varbio_male_30to35) sum_varbio_male_30to35,
sum(varbio_male_35to40) sum_varbio_male_35to40,
sum(varbio_male_40to50) sum_varbio_male_40to50,
sum(varbio_male_40to45) sum_varbio_male_40to45,
sum(varbio_male_45to50) sum_varbio_male_45to50,
sum(varbio_male_50to60) sum_varbio_male_50to60,
sum(varbio_male_30to50) sum_varbio_male_30to50,
sum(varbio_male_30to60) sum_varbio_male_30to60,
sum(varbio_male_le103) sum_varbio_male_le103,
sum(varbio_female_le30) sum_varbio_female_le30,
sum(varbio_female_le40) sum_varbio_female_le40,
sum(varbio_female_le50) sum_varbio_female_le50,
sum(varbio_female_le60) sum_varbio_female_le60,
sum(varbio_female_30to40) sum_varbio_female_30to40,
sum(varbio_female_30to35) sum_varbio_female_30to35,
sum(varbio_female_35to40) sum_varbio_female_35to40,
sum(varbio_female_40to50) sum_varbio_female_40to50,
sum(varbio_female_40to45) sum_varbio_female_40to45,
sum(varbio_female_45to50) sum_varbio_female_45to50,
sum(varbio_female_50to60) sum_varbio_female_50to60,
sum(varbio_female_30to50) sum_varbio_female_30to50,
sum(varbio_female_30to60) sum_varbio_female_30to60,
sum(varbio_female_immature) sum_varbio_female_immature
from cb_variancebio_sizegroup
group by survey_year
order by survey_year;

drop table cb_pop_sizegroup_cv;

create table cb_pop_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_pop_male_le30 <> 0
	 then ((sqrt(sum_varpop_male_le30))/sum_pop_male_le30)
	 else 0
	 end) cv_pop_male_le30,
(CASE
	 when sum_pop_male_le40 <> 0
	 then ((sqrt(sum_varpop_male_le40))/sum_pop_male_le40)
	 else 0
	 end) cv_pop_male_le40,   
(CASE
	 when sum_pop_male_le50 <> 0
	 then ((sqrt(sum_varpop_male_le50))/sum_pop_male_le50)
	 else 0
	 end) cv_pop_male_le50,
(CASE
	 when sum_pop_male_le60 <> 0
	 then ((sqrt(sum_varpop_male_le60))/sum_pop_male_le60)
	 else 0
	 end) cv_pop_male_le60,
(CASE
	 when sum_pop_male_30to40 <> 0
	 then ((sqrt(sum_varpop_male_30to40))/sum_pop_male_30to40)
	 else 0
	 end) cv_pop_male_30to40,
(CASE
	 when sum_pop_male_30to35 <> 0
	 then ((sqrt(sum_varpop_male_30to35))/sum_pop_male_30to35)
	 else 0
	 end) cv_pop_male_30to35,
(CASE
	 when sum_pop_male_35to40 <> 0
	 then ((sqrt(sum_varpop_male_35to40))/sum_pop_male_35to40)
	 else 0
	 end) cv_pop_male_35to40,
(CASE
	 when sum_pop_male_40to50 <> 0
	 then ((sqrt(sum_varpop_male_40to50))/sum_pop_male_40to50)
	 else 0
	 end) cv_pop_male_40to50,
(CASE
	 when sum_pop_male_40to45 <> 0
	 then ((sqrt(sum_varpop_male_40to45))/sum_pop_male_40to45)
	 else 0
	 end) cv_pop_male_40to45,
(CASE
	 when sum_pop_male_45to50 <> 0
	 then ((sqrt(sum_varpop_male_45to50))/sum_pop_male_45to50)
	 else 0
	 end) cv_pop_male_45to50,
(CASE
	 when sum_pop_male_50to60 <> 0
	 then ((sqrt(sum_varpop_male_50to60))/sum_pop_male_50to60)
	 else 0
	 end) cv_pop_male_50to60,
(CASE
	 when sum_pop_male_30to50 <> 0
	 then ((sqrt(sum_varpop_male_30to50))/sum_pop_male_30to50)
	 else 0
	 end) cv_pop_male_30to50,   
(CASE
	 when sum_pop_male_30to60 <> 0
	 then ((sqrt(sum_varpop_male_30to60))/sum_pop_male_30to60)
	 else 0
	 end) cv_pop_male_30to60,
(CASE
	 when sum_pop_male_le103 <> 0
	 then ((sqrt(sum_varpop_male_le103))/sum_pop_male_le103)
	 else 0
	 end) cv_pop_male_le103,
(CASE
	 when sum_pop_female_le30 <> 0
	 then ((sqrt(sum_varpop_female_le30))/sum_pop_female_le30)
	 else 0
	 end) cv_pop_female_le30,
(CASE
	 when sum_pop_female_le40 <> 0
	 then ((sqrt(sum_varpop_female_le40))/sum_pop_female_le40)
	 else 0
	 end) cv_pop_female_le40,   
(CASE
	 when sum_pop_female_le50 <> 0
	 then ((sqrt(sum_varpop_female_le50))/sum_pop_female_le50)
	 else 0
	 end) cv_pop_female_le50,
(CASE
	 when sum_pop_female_le60 <> 0
	 then ((sqrt(sum_varpop_female_le60))/sum_pop_female_le60)
	 else 0
	 end) cv_pop_female_le60,
(CASE
	 when sum_pop_female_30to40 <> 0
	 then ((sqrt(sum_varpop_female_30to40))/sum_pop_female_30to40)
	 else 0
	 end) cv_pop_female_30to40,
(CASE
	 when sum_pop_female_30to35 <> 0
	 then ((sqrt(sum_varpop_female_30to35))/sum_pop_female_30to35)
	 else 0
	 end) cv_pop_female_30to35,
(CASE
	 when sum_pop_female_35to40 <> 0
	 then ((sqrt(sum_varpop_female_35to40))/sum_pop_female_35to40)
	 else 0
	 end) cv_pop_female_35to40,
(CASE
	 when sum_pop_female_40to50 <> 0
	 then ((sqrt(sum_varpop_female_40to50))/sum_pop_female_40to50)
	 else 0
	 end) cv_pop_female_40to50,
(CASE
	 when sum_pop_female_40to45 <> 0
	 then ((sqrt(sum_varpop_female_40to45))/sum_pop_female_40to45)
	 else 0
	 end) cv_pop_female_40to45,
(CASE
	 when sum_pop_female_45to50 <> 0
	 then ((sqrt(sum_varpop_female_45to50))/sum_pop_female_45to50)
	 else 0
	 end) cv_pop_female_45to50,
(CASE
	 when sum_pop_female_50to60 <> 0
	 then ((sqrt(sum_varpop_female_50to60))/sum_pop_female_50to60)
	 else 0
	 end) cv_pop_female_50to60,
(CASE
	 when sum_pop_female_30to50 <> 0
	 then ((sqrt(sum_varpop_female_30to50))/sum_pop_female_30to50)
	 else 0
	 end) cv_pop_female_30to50,   
(CASE
	 when sum_pop_female_30to60 <> 0
	 then ((sqrt(sum_varpop_female_30to60))/sum_pop_female_30to60)
	 else 0
	 end) cv_pop_female_30to60,
(CASE
	 when sum_pop_female_immature <> 0
	 then ((sqrt(sum_varpop_female_immature))/sum_pop_female_immature)
	 else 0
	 end) cv_pop_female_immature 	 	 	 	  

from cb_varpop_sizegroup_sum a, cb_popall_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;	 	 	 

drop table cb_bio_sizegroup_cv;

create table cb_bio_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_bio_male_le30 <> 0
	 then ((sqrt(sum_varbio_male_le30))/sum_bio_male_le30)
	 else 0
	 end) cv_bio_male_le30,
(CASE
	 when sum_bio_male_le40 <> 0
	 then ((sqrt(sum_varbio_male_le40))/sum_bio_male_le40)
	 else 0
	 end) cv_bio_male_le40,   
(CASE
	 when sum_bio_male_le50 <> 0
	 then ((sqrt(sum_varbio_male_le50))/sum_bio_male_le50)
	 else 0
	 end) cv_bio_male_le50,
(CASE
	 when sum_bio_male_le60 <> 0
	 then ((sqrt(sum_varbio_male_le60))/sum_bio_male_le60)
	 else 0
	 end) cv_bio_male_le60,
(CASE
	 when sum_bio_male_30to40 <> 0
	 then ((sqrt(sum_varbio_male_30to40))/sum_bio_male_30to40)
	 else 0
	 end) cv_bio_male_30to40,
(CASE
	 when sum_bio_male_30to35 <> 0
	 then ((sqrt(sum_varbio_male_30to35))/sum_bio_male_30to35)
	 else 0
	 end) cv_bio_male_30to35,
(CASE
	 when sum_bio_male_35to40 <> 0
	 then ((sqrt(sum_varbio_male_35to40))/sum_bio_male_35to40)
	 else 0
	 end) cv_bio_male_35to40,
(CASE
	 when sum_bio_male_40to50 <> 0
	 then ((sqrt(sum_varbio_male_40to50))/sum_bio_male_40to50)
	 else 0
	 end) cv_bio_male_40to50,
(CASE
	 when sum_bio_male_40to45 <> 0
	 then ((sqrt(sum_varbio_male_40to45))/sum_bio_male_40to45)
	 else 0
	 end) cv_bio_male_40to45,
(CASE
	 when sum_bio_male_45to50 <> 0
	 then ((sqrt(sum_varbio_male_45to50))/sum_bio_male_45to50)
	 else 0
	 end) cv_bio_male_45to50,
(CASE
	 when sum_bio_male_50to60 <> 0
	 then ((sqrt(sum_varbio_male_50to60))/sum_bio_male_50to60)
	 else 0
	 end) cv_bio_male_50to60,
(CASE
	 when sum_bio_male_30to50 <> 0
	 then ((sqrt(sum_varbio_male_30to50))/sum_bio_male_30to50)
	 else 0
	 end) cv_bio_male_30to50,
(CASE
	 when sum_bio_male_30to60 <> 0
	 then ((sqrt(sum_varbio_male_30to60))/sum_bio_male_30to60)
	 else 0
	 end) cv_bio_male_30to60,   
(CASE
	 when sum_bio_male_le103 <> 0
	 then ((sqrt(sum_varbio_male_le103))/sum_bio_male_le103)
	 else 0
	 end) cv_bio_male_le103,
(CASE
	 when sum_bio_female_le30 <> 0
	 then ((sqrt(sum_varbio_female_le30))/sum_bio_female_le30)
	 else 0
	 end) cv_bio_female_le30,
(CASE
	 when sum_bio_female_le40 <> 0
	 then ((sqrt(sum_varbio_female_le40))/sum_bio_female_le40)
	 else 0
	 end) cv_bio_female_le40,   
(CASE
	 when sum_bio_female_le50 <> 0
	 then ((sqrt(sum_varbio_female_le50))/sum_bio_female_le50)
	 else 0
	 end) cv_bio_female_le50,
(CASE
	 when sum_bio_female_le60 <> 0
	 then ((sqrt(sum_varbio_female_le60))/sum_bio_female_le60)
	 else 0
	 end) cv_bio_female_le60,
(CASE
	 when sum_bio_female_30to40 <> 0
	 then ((sqrt(sum_varbio_female_30to40))/sum_bio_female_30to40)
	 else 0
	 end) cv_bio_female_30to40,
(CASE
	 when sum_bio_female_30to35 <> 0
	 then ((sqrt(sum_varbio_female_30to35))/sum_bio_female_30to35)
	 else 0
	 end) cv_bio_female_30to35,
(CASE
	 when sum_bio_female_35to40 <> 0
	 then ((sqrt(sum_varbio_female_35to40))/sum_bio_female_35to40)
	 else 0
	 end) cv_bio_female_35to40,
(CASE
	 when sum_bio_female_40to50 <> 0
	 then ((sqrt(sum_varbio_female_40to50))/sum_bio_female_40to50)
	 else 0
	 end) cv_bio_female_40to50,
(CASE
	 when sum_bio_female_40to45 <> 0
	 then ((sqrt(sum_varbio_female_40to45))/sum_bio_female_40to45)
	 else 0
	 end) cv_bio_female_40to45,
(CASE
	 when sum_bio_female_45to50 <> 0
	 then ((sqrt(sum_varbio_female_45to50))/sum_bio_female_45to50)
	 else 0
	 end) cv_bio_female_45to50,
(CASE
	 when sum_bio_female_50to60 <> 0
	 then ((sqrt(sum_varbio_female_50to60))/sum_bio_female_50to60)
	 else 0
	 end) cv_bio_female_50to60,
(CASE
	 when sum_bio_female_30to50 <> 0
	 then ((sqrt(sum_varbio_female_30to50))/sum_bio_female_30to50)
	 else 0
	 end) cv_bio_female_30to50,
(CASE
	 when sum_bio_female_30to60 <> 0
	 then ((sqrt(sum_varbio_female_30to60))/sum_bio_female_30to60)
	 else 0
	 end) cv_bio_female_30to60,

(CASE
	 when sum_bio_female_immature <> 0
	 then ((sqrt(sum_varbio_female_immature))/sum_bio_female_immature)
	 else 0
	 end) cv_bio_female_immature 	 	 	 	  
from cb_varbio_sizegroup_sum a, cb_bioall_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;


-- CI calcs

create or replace view cb_sizegroup_stderr_pop as
select distinct survey_year,
(sqrt(sum_varpop_male_le30)) stderr_pop_male_le30,
(sqrt(sum_varpop_male_le40)) stderr_pop_male_le40,
(sqrt(sum_varpop_male_le50)) stderr_pop_male_le50,
(sqrt(sum_varpop_male_le60)) stderr_pop_male_le60,
(sqrt(sum_varpop_male_30to40)) stderr_pop_male_30to40,
(sqrt(sum_varpop_male_30to35)) stderr_pop_male_30to35,
(sqrt(sum_varpop_male_35to40)) stderr_pop_male_35to40,
(sqrt(sum_varpop_male_40to50)) stderr_pop_male_40to50,
(sqrt(sum_varpop_male_40to45)) stderr_pop_male_40to45,
(sqrt(sum_varpop_male_45to50)) stderr_pop_male_45to50,
(sqrt(sum_varpop_male_50to60)) stderr_pop_male_50to60,
(sqrt(sum_varpop_male_30to50)) stderr_pop_male_30to50,
(sqrt(sum_varpop_male_30to60)) stderr_pop_male_30to60,
(sqrt(sum_varpop_male_le103)) stderr_pop_male_le103,
(sqrt(sum_varpop_female_le30)) stderr_pop_female_le30,
(sqrt(sum_varpop_female_le40)) stderr_pop_female_le40,
(sqrt(sum_varpop_female_le50)) stderr_pop_female_le50,
(sqrt(sum_varpop_female_le60)) stderr_pop_female_le60,
(sqrt(sum_varpop_female_30to40)) stderr_pop_female_30to40,
(sqrt(sum_varpop_female_30to35)) stderr_pop_female_30to35,
(sqrt(sum_varpop_female_35to40)) stderr_pop_female_35to40,
(sqrt(sum_varpop_female_40to50)) stderr_pop_female_40to50,
(sqrt(sum_varpop_female_40to45)) stderr_pop_female_40to45,
(sqrt(sum_varpop_female_45to50)) stderr_pop_female_45to50,
(sqrt(sum_varpop_female_50to60)) stderr_pop_female_50to60,
(sqrt(sum_varpop_female_30to50)) stderr_pop_female_30to50,
(sqrt(sum_varpop_female_30to60)) stderr_pop_female_30to60,
(sqrt(sum_varpop_female_immature)) stderr_pop_female_immature
from cb_varpop_sizegroup_sum;

create or replace view cb_sizegroup_stderr_bio as
select distinct survey_year,
(sqrt(sum_varbio_male_le30)) stderr_bio_male_le30,
(sqrt(sum_varbio_male_le40)) stderr_bio_male_le40,
(sqrt(sum_varbio_male_le50)) stderr_bio_male_le50,
(sqrt(sum_varbio_male_le60)) stderr_bio_male_le60,
(sqrt(sum_varbio_male_30to40)) stderr_bio_male_30to40,
(sqrt(sum_varbio_male_30to35)) stderr_bio_male_30to35,
(sqrt(sum_varbio_male_35to40)) stderr_bio_male_35to40,
(sqrt(sum_varbio_male_40to50)) stderr_bio_male_40to50,
(sqrt(sum_varbio_male_40to45)) stderr_bio_male_40to45,
(sqrt(sum_varbio_male_45to50)) stderr_bio_male_45to50,
(sqrt(sum_varbio_male_50to60)) stderr_bio_male_50to60,
(sqrt(sum_varbio_male_30to50)) stderr_bio_male_30to50,
(sqrt(sum_varbio_male_30to60)) stderr_bio_male_30to60,
(sqrt(sum_varbio_male_le103)) stderr_bio_male_le103,
(sqrt(sum_varbio_female_le30)) stderr_bio_female_le30,
(sqrt(sum_varbio_female_le40)) stderr_bio_female_le40,
(sqrt(sum_varbio_female_le50)) stderr_bio_female_le50,
(sqrt(sum_varbio_female_le60)) stderr_bio_female_le60,
(sqrt(sum_varbio_female_30to40)) stderr_bio_female_30to40,
(sqrt(sum_varbio_female_30to35)) stderr_bio_female_30to35,
(sqrt(sum_varbio_female_35to40)) stderr_bio_female_35to40,
(sqrt(sum_varbio_female_40to50)) stderr_bio_female_40to50,
(sqrt(sum_varbio_female_40to45)) stderr_bio_female_40to45,
(sqrt(sum_varbio_female_45to50)) stderr_bio_female_45to50,
(sqrt(sum_varbio_female_50to60)) stderr_bio_female_50to60,
(sqrt(sum_varbio_female_30to50)) stderr_bio_female_30to50,
(sqrt(sum_varbio_female_30to60)) stderr_bio_female_30to60,
(sqrt(sum_varbio_female_immature)) stderr_bio_female_immature
from cb_varbio_sizegroup_sum;

drop table cb_sizegroup_confidence_pop;

create table cb_sizegroup_confidence_pop as
select distinct survey_year,
(1.96 * stderr_pop_male_le30) ci_pop_male_le30,
(1.96 * stderr_pop_male_le40) ci_pop_male_le40,
(1.96 * stderr_pop_male_le50) ci_pop_male_le50,
(1.96 * stderr_pop_male_le60) ci_pop_male_le60,
(1.96 * stderr_pop_male_30to40) ci_pop_male_30to40,
(1.96 * stderr_pop_male_30to35) ci_pop_male_30to35,
(1.96 * stderr_pop_male_35to40) ci_pop_male_35to40,
(1.96 * stderr_pop_male_40to50) ci_pop_male_40to50,
(1.96 * stderr_pop_male_40to45) ci_pop_male_40to45,
(1.96 * stderr_pop_male_45to50) ci_pop_male_45to50,
(1.96 * stderr_pop_male_50to60) ci_pop_male_50to60,
(1.96 * stderr_pop_male_30to50) ci_pop_male_30to50,
(1.96 * stderr_pop_male_30to60) ci_pop_male_30to60,
(1.96 * stderr_pop_male_le103) ci_pop_male_le103,
(1.96 * stderr_pop_female_le30) ci_pop_female_le30,
(1.96 * stderr_pop_female_le40) ci_pop_female_le40,
(1.96 * stderr_pop_female_le50) ci_pop_female_le50,
(1.96 * stderr_pop_female_le60) ci_pop_female_le60,
(1.96 * stderr_pop_female_30to40) ci_pop_female_30to40,
(1.96 * stderr_pop_female_30to35) ci_pop_female_30to35,
(1.96 * stderr_pop_female_35to40) ci_pop_female_35to40,
(1.96 * stderr_pop_female_40to50) ci_pop_female_40to50,
(1.96 * stderr_pop_female_40to45) ci_pop_female_40to45,
(1.96 * stderr_pop_female_45to50) ci_pop_female_45to50,
(1.96 * stderr_pop_female_50to60) ci_pop_female_50to60,
(1.96 * stderr_pop_female_30to50) ci_pop_female_30to50,
(1.96 * stderr_pop_female_30to60) ci_pop_female_30to60,
(1.96 * stderr_pop_female_immature) ci_pop_female_immature
from cb_sizegroup_stderr_pop;

drop table cb_sizegroup_confidence_bio;

create table cb_sizegroup_confidence_bio as
select distinct survey_year,
((1.96 * stderr_bio_male_le30)) ci_bio_male_le30,
((1.96 * stderr_bio_male_le40)) ci_bio_male_le40,
((1.96 * stderr_bio_male_le50)) ci_bio_male_le50,
((1.96 * stderr_bio_male_le60)) ci_bio_male_le60,
((1.96 * stderr_bio_male_30to40)) ci_bio_male_30to40,
((1.96 * stderr_bio_male_30to35)) ci_bio_male_30to35,
((1.96 * stderr_bio_male_35to40)) ci_bio_male_35to40,
((1.96 * stderr_bio_male_40to50)) ci_bio_male_40to50,
((1.96 * stderr_bio_male_40to45)) ci_bio_male_40to45,
((1.96 * stderr_bio_male_45to50)) ci_bio_male_45to50,
((1.96 * stderr_bio_male_50to60)) ci_bio_male_50to60,
((1.96 * stderr_bio_male_30to50)) ci_bio_male_30to50,
((1.96 * stderr_bio_male_30to60)) ci_bio_male_30to60,
((1.96 * stderr_bio_male_le103)) ci_bio_male_le103,
((1.96 * stderr_bio_female_le30)) ci_bio_female_le30,
((1.96 * stderr_bio_female_le40)) ci_bio_female_le40,
((1.96 * stderr_bio_female_le50)) ci_bio_female_le50,
((1.96 * stderr_bio_female_le60)) ci_bio_female_le60,
((1.96 * stderr_bio_female_30to40)) ci_bio_female_30to40,
((1.96 * stderr_bio_female_30to35)) ci_bio_female_30to35,
((1.96 * stderr_bio_female_35to40)) ci_bio_female_35to40,
((1.96 * stderr_bio_female_40to50)) ci_bio_female_40to50,
((1.96 * stderr_bio_female_40to45)) ci_bio_female_40to45,
((1.96 * stderr_bio_female_45to50)) ci_bio_female_45to50,
((1.96 * stderr_bio_female_50to60)) ci_bio_female_50to60,
((1.96 * stderr_bio_female_30to50)) ci_bio_female_30to50,
((1.96 * stderr_bio_female_30to60)) ci_bio_female_30to60,
((1.96 * stderr_bio_female_immature)) ci_bio_female_immature
from cb_sizegroup_stderr_bio;


-- Final output for stocks

-- All Districts Combined

drop table cb_ebs_bio_juvenile;

create table cb_ebs_bio_juvenile as
select a.survey_year,
sum_bio_male_le30 biomass_male_le30,cv_bio_male_le30 cv_biomass_male_le30,ci_bio_male_le30 ci_biomass_male_le30,
sum_bio_male_le40 biomass_male_le40,cv_bio_male_le40 cv_biomass_male_le40,ci_bio_male_le40 ci_biomass_male_le40,
sum_bio_male_le50 biomass_male_le50,cv_bio_male_le50 cv_biomass_male_le50,ci_bio_male_le50 ci_biomass_male_le50,
sum_bio_male_le60 biomass_male_le60,cv_bio_male_le60 cv_biomass_male_le60,ci_bio_male_le60 ci_biomass_male_le60,
sum_bio_male_30to35 biomass_male_30to35,cv_bio_male_30to35 cv_biomass_male_30to35,ci_bio_male_30to35 ci_biomass_male_30to35,
sum_bio_male_35to40 biomass_male_35to40,cv_bio_male_35to40 cv_biomass_male_35to40,ci_bio_male_35to40 ci_biomass_male_35to40,
sum_bio_male_30to40 biomass_male_30to40,cv_bio_male_30to40 cv_biomass_male_30to40,ci_bio_male_30to40 ci_biomass_male_30to40,
sum_bio_male_40to45 biomass_male_40to45,cv_bio_male_40to45 cv_biomass_male_40to45,ci_bio_male_40to45 ci_biomass_male_40to45,
sum_bio_male_45to50 biomass_male_45to50,cv_bio_male_45to50 cv_biomass_male_45to50,ci_bio_male_45to50 ci_biomass_male_45to50,
sum_bio_male_40to50 biomass_male_40to50,cv_bio_male_40to50 cv_biomass_male_40to50,ci_bio_male_40to50 ci_biomass_male_40to50,
sum_bio_male_50to60 biomass_male_50to60,cv_bio_male_50to60 cv_biomass_male_50to60,ci_bio_male_50to60 ci_biomass_male_50to60,
sum_bio_male_30to50 biomass_male_30to50,cv_bio_male_30to50 cv_biomass_male_30to50,ci_bio_male_30to50 ci_biomass_male_30to50,
sum_bio_male_30to60 biomass_male_30to60,cv_bio_male_30to60 cv_biomass_male_30to60,ci_bio_male_30to60 ci_biomass_male_30to60,
sum_bio_male_le103 biomass_male_le103,cv_bio_male_le103 cv_biomass_male_le103,ci_bio_male_le103 ci_biomass_male_le103,
sum_bio_female_le30 biomass_female_le30,cv_bio_female_le30 cv_biomass_female_le30,ci_bio_female_le30 ci_biomass_female_le30,
sum_bio_female_le40 biomass_female_le40,cv_bio_female_le40 cv_biomass_female_le40,ci_bio_female_le40 ci_biomass_female_le40,
sum_bio_female_le50 biomass_female_le50,cv_bio_female_le50 cv_biomass_female_le50,ci_bio_female_le50 ci_biomass_female_le50,
sum_bio_female_le60 biomass_female_le60,cv_bio_female_le60 cv_biomass_female_le60,ci_bio_female_le60 ci_biomass_female_le60,
sum_bio_female_30to35 biomass_female_30to35,cv_bio_female_30to35 cv_biomass_female_30to35,ci_bio_female_30to35 ci_biomass_female_30to35,
sum_bio_female_35to40 biomass_female_35to40,cv_bio_female_35to40 cv_biomass_female_35to40,ci_bio_female_35to40 ci_biomass_female_35to40,
sum_bio_female_30to40 biomass_female_30to40,cv_bio_female_30to40 cv_biomass_female_30to40,ci_bio_female_30to40 ci_biomass_female_30to40,
sum_bio_female_40to45 biomass_female_40to45,cv_bio_female_40to45 cv_biomass_female_40to45,ci_bio_female_40to45 ci_biomass_female_40to45,
sum_bio_female_45to50 biomass_female_45to50,cv_bio_female_45to50 cv_biomass_female_45to50,ci_bio_female_45to50 ci_biomass_female_45to50,
sum_bio_female_40to50 biomass_female_40to50,cv_bio_female_40to50 cv_biomass_female_40to50,ci_bio_female_40to50 ci_biomass_female_40to50,
sum_bio_female_50to60 biomass_female_50to60,cv_bio_female_50to60 cv_biomass_female_50to60,ci_bio_female_50to60 ci_biomass_female_50to60,
sum_bio_female_30to50 biomass_female_30to50,cv_bio_female_30to50 cv_biomass_female_30to50,ci_bio_female_30to50 ci_biomass_female_30to50,
sum_bio_female_30to60 biomass_female_30to60,cv_bio_female_30to60 cv_biomass_female_30to60,ci_bio_female_30to60 ci_biomass_female_30to60,
sum_bio_female_immature biomass_female_immature,cv_bio_female_immature cv_biomass_female_immature,ci_bio_female_immature ci_biomass_female_immature
from cb_bioall_sizegroup a, cb_sizegroup_confidence_bio b,cb_bio_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;


drop table cb_ebs_pop_juvenile;

create table cb_ebs_pop_juvenile as
select a.survey_year,
sum_pop_male_le30 num_male_le30,cv_pop_male_le30 cv_num_male_le30,ci_pop_male_le30 ci_num_male_le30,
sum_pop_male_le40 num_male_le40,cv_pop_male_le40 cv_num_male_le40,ci_pop_male_le40 ci_num_male_le40,
sum_pop_male_le50 num_male_le50,cv_pop_male_le50 cv_num_male_le50,ci_pop_male_le50 ci_num_male_le50,
sum_pop_male_le60 num_male_le60,cv_pop_male_le60 cv_num_male_le60,ci_pop_male_le60 ci_num_male_le60,

sum_pop_male_30to35 num_male_30to35,cv_pop_male_30to35 cv_num_male_30to35,ci_pop_male_30to35 ci_num_male_30to35,
sum_pop_male_35to40 num_male_35to40,cv_pop_male_35to40 cv_num_male_35to40,ci_pop_male_35to40 ci_num_male_35to40,

sum_pop_male_30to40 num_male_30to40,cv_pop_male_30to40 cv_num_male_30to40,ci_pop_male_30to40 ci_num_male_30to40,

sum_pop_male_40to45 num_male_40to45,cv_pop_male_40to45 cv_num_male_40to45,ci_pop_male_40to45 ci_num_male_40to45,
sum_pop_male_45to50 num_male_45to50,cv_pop_male_45to50 cv_num_male_45to50,ci_pop_male_45to50 ci_num_male_45to50,

sum_pop_male_40to50 num_male_40to50,cv_pop_male_40to50 cv_num_male_40to50,ci_pop_male_40to50 ci_num_male_40to50,
sum_pop_male_50to60 num_male_50to60,cv_pop_male_50to60 cv_num_male_50to60,ci_pop_male_50to60 ci_num_male_50to60,
sum_pop_male_30to50 num_male_30to50,cv_pop_male_30to50 cv_num_male_30to50,ci_pop_male_30to50 ci_num_male_30to50,
sum_pop_male_30to60 num_male_30to60,cv_pop_male_30to60 cv_num_male_30to60,ci_pop_male_30to60 ci_num_male_30to60,
sum_pop_male_le103 num_male_le103,cv_pop_male_le103 cv_num_male_le103,ci_pop_male_le103 ci_num_male_le103,
sum_pop_female_le30 num_female_le30,cv_pop_female_le30 cv_num_female_le30,ci_pop_female_le30 ci_num_female_le30,
sum_pop_female_le40 num_female_le40,cv_pop_female_le40 cv_num_female_le40,ci_pop_female_le40 ci_num_female_le40,
sum_pop_female_le50 num_female_le50,cv_pop_female_le50 cv_num_female_le50,ci_pop_female_le50 ci_num_female_le50,
sum_pop_female_le60 num_female_le60,cv_pop_female_le60 cv_num_female_le60,ci_pop_female_le60 ci_num_female_le60,

sum_pop_female_30to35 num_female_30to35,cv_pop_female_30to35 cv_num_female_30to35,ci_pop_female_30to35 ci_num_female_30to35,
sum_pop_female_35to40 num_female_35to40,cv_pop_female_35to40 cv_num_female_35to40,ci_pop_female_35to40 ci_num_female_35to40,

sum_pop_female_30to40 num_female_30to40,cv_pop_female_30to40 cv_num_female_30to40,ci_pop_female_30to40 ci_num_female_30to40,

sum_pop_female_40to45 num_female_40to45,cv_pop_female_40to45 cv_num_female_40to45,ci_pop_female_40to45 ci_num_female_40to45,
sum_pop_female_45to50 num_female_45to50,cv_pop_female_45to50 cv_num_female_45to50,ci_pop_female_45to50 ci_num_female_45to50,

sum_pop_female_40to50 num_female_40to50,cv_pop_female_40to50 cv_num_female_40to50,ci_pop_female_40to50 ci_num_female_40to50,
sum_pop_female_50to60 num_female_50to60,cv_pop_female_50to60 cv_num_female_50to60,ci_pop_female_50to60 ci_num_female_50to60,
sum_pop_female_30to50 num_female_30to50,cv_pop_female_30to50 cv_num_female_30to50,ci_pop_female_30to50 ci_num_female_30to50,
sum_pop_female_30to60 num_female_30to60,cv_pop_female_30to60 cv_num_female_30to60,ci_pop_female_30to60 ci_num_female_30to60,
sum_pop_female_immature num_female_immature,cv_pop_female_immature cv_num_female_immature,ci_pop_female_immature ci_num_female_immature
from cb_popall_sizegroup a, cb_sizegroup_confidence_pop b,cb_pop_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
order by a.survey_year;


-- FOR QUADRANTS DELETE/SHUTDOWN East-166 and West-166
--  east of 166
drop table e166cb_pop_sizegroup;

create table e166cb_pop_sizegroup as
select survey_year,
sum(pop_male_le30) sum_pop_male_le30,
sum(pop_male_le40) sum_pop_male_le40,
sum(pop_male_le50) sum_pop_male_le50,
sum(pop_male_le60) sum_pop_male_le60,
sum(pop_male_30to35) sum_pop_male_30to35,
sum(pop_male_35to40) sum_pop_male_35to40,
sum(pop_male_30to40) sum_pop_male_30to40,
sum(pop_male_40to45) sum_pop_male_40to45,
sum(pop_male_45to50) sum_pop_male_45to50,
sum(pop_male_40to50) sum_pop_male_40to50,
sum(pop_male_50to60) sum_pop_male_50to60,
sum(pop_male_30to50) sum_pop_male_30to50,
sum(pop_male_30to60) sum_pop_male_30to60,
sum(pop_male_le103) sum_pop_male_le103,
sum(pop_female_le30) sum_pop_female_le30,
sum(pop_female_le40) sum_pop_female_le40,
sum(pop_female_le50) sum_pop_female_le50,
sum(pop_female_le60) sum_pop_female_le60,
sum(pop_female_30to35) sum_pop_female_30to35,
sum(pop_female_35to40) sum_pop_female_35to40,
sum(pop_female_30to40) sum_pop_female_30to40,
sum(pop_female_40to45) sum_pop_female_40to45,
sum(pop_female_45to50) sum_pop_female_45to50,
sum(pop_female_40to50) sum_pop_female_40to50,
sum(pop_female_50to60) sum_pop_female_50to60,
sum(pop_female_30to50) sum_pop_female_30to50,
sum(pop_female_30to60) sum_pop_female_30to60,
sum(pop_female_immature) sum_pop_female_immature
from cb_popbystratum_sizegroup
where district = 'East 166'
group by survey_year
order by survey_year;

drop table e166cb_bio_sizegroup;

create table e166cb_bio_sizegroup as
select survey_year,
sum(bio_male_le30) sum_bio_male_le30,
sum(bio_male_le40) sum_bio_male_le40,
sum(bio_male_le50) sum_bio_male_le50,
sum(bio_male_le60) sum_bio_male_le60,
sum(bio_male_30to35) sum_bio_male_30to35,
sum(bio_male_35to40) sum_bio_male_35to40,
sum(bio_male_30to40) sum_bio_male_30to40,
sum(bio_male_40to45) sum_bio_male_40to45,
sum(bio_male_45to50) sum_bio_male_45to50,
sum(bio_male_40to50) sum_bio_male_40to50,
sum(bio_male_50to60) sum_bio_male_50to60,
sum(bio_male_30to50) sum_bio_male_30to50,
sum(bio_male_30to60) sum_bio_male_30to60,
sum(bio_male_le103) sum_bio_male_le103,
sum(bio_female_le30) sum_bio_female_le30,
sum(bio_female_le40) sum_bio_female_le40,
sum(bio_female_le50) sum_bio_female_le50,
sum(bio_female_le60) sum_bio_female_le60,
sum(bio_female_30to35) sum_bio_female_30to35,
sum(bio_female_35to40) sum_bio_female_35to40,
sum(bio_female_30to40) sum_bio_female_30to40,
sum(bio_female_40to45) sum_bio_female_40to45,
sum(bio_female_45to50) sum_bio_female_45to50,
sum(bio_female_40to50) sum_bio_female_40to50,
sum(bio_female_50to60) sum_bio_female_50to60,
sum(bio_female_30to50) sum_bio_female_30to50,
sum(bio_female_30to60) sum_bio_female_30to60,
sum(bio_female_immature) sum_bio_female_immature
from cb_biobystratum_sizegroup
where district = 'East 166'
group by survey_year
order by survey_year;



drop table e166cb_varpop_sizegroup_sum;

create table e166cb_varpop_sizegroup_sum as
select distinct survey_year,
sum(varpop_male_le30) sum_varpop_male_le30,
sum(varpop_male_le40) sum_varpop_male_le40,
sum(varpop_male_le50) sum_varpop_male_le50,
sum(varpop_male_le60) sum_varpop_male_le60,
sum(varpop_male_30to35) sum_varpop_male_30to35,
sum(varpop_male_35to40) sum_varpop_male_35to40,
sum(varpop_male_30to40) sum_varpop_male_30to40,
sum(varpop_male_40to45) sum_varpop_male_40to45,
sum(varpop_male_45to50) sum_varpop_male_45to50,
sum(varpop_male_40to50) sum_varpop_male_40to50,
sum(varpop_male_50to60) sum_varpop_male_50to60,
sum(varpop_male_30to50) sum_varpop_male_30to50,
sum(varpop_male_30to60) sum_varpop_male_30to60,
sum(varpop_male_le103) sum_varpop_male_le103,
sum(varpop_female_le30) sum_varpop_female_le30,
sum(varpop_female_le40) sum_varpop_female_le40,
sum(varpop_female_le50) sum_varpop_female_le50,
sum(varpop_female_le60) sum_varpop_female_le60,
sum(varpop_female_30to35) sum_varpop_female_30to35,
sum(varpop_female_35to40) sum_varpop_female_35to40,
sum(varpop_female_30to40) sum_varpop_female_30to40,
sum(varpop_female_40to45) sum_varpop_female_40to45,
sum(varpop_female_45to50) sum_varpop_female_45to50,
sum(varpop_female_40to50) sum_varpop_female_40to50,
sum(varpop_female_50to60) sum_varpop_female_50to60,
sum(varpop_female_30to50) sum_varpop_female_30to50,
sum(varpop_female_30to60) sum_varpop_female_30to60,
sum(varpop_female_immature) sum_varpop_female_immature
from cb_variancepop_sizegroup
where district = 'East 166'
group by survey_year
order by survey_year;

drop table e166cb_varbio_sizegroup_sum;

create table e166cb_varbio_sizegroup_sum as
select distinct survey_year,
sum(varbio_male_le30) sum_varbio_male_le30,
sum(varbio_male_le40) sum_varbio_male_le40,
sum(varbio_male_le50) sum_varbio_male_le50,
sum(varbio_male_le60) sum_varbio_male_le60,
sum(varbio_male_30to35) sum_varbio_male_30to35,
sum(varbio_male_35to40) sum_varbio_male_35to40,
sum(varbio_male_30to40) sum_varbio_male_30to40,
sum(varbio_male_40to45) sum_varbio_male_40to45,
sum(varbio_male_45to50) sum_varbio_male_45to50,
sum(varbio_male_40to50) sum_varbio_male_40to50,
sum(varbio_male_50to60) sum_varbio_male_50to60,
sum(varbio_male_30to50) sum_varbio_male_30to50,
sum(varbio_male_30to60) sum_varbio_male_30to60,
sum(varbio_male_le103) sum_varbio_male_le103,
sum(varbio_female_le30) sum_varbio_female_le30,
sum(varbio_female_le40) sum_varbio_female_le40,
sum(varbio_female_le50) sum_varbio_female_le50,
sum(varbio_female_le60) sum_varbio_female_le60,
sum(varbio_female_30to35) sum_varbio_female_30to35,
sum(varbio_female_35to40) sum_varbio_female_35to40,
sum(varbio_female_30to40) sum_varbio_female_30to40,
sum(varbio_female_40to45) sum_varbio_female_40to45,
sum(varbio_female_45to50) sum_varbio_female_45to50,
sum(varbio_female_40to50) sum_varbio_female_40to50,
sum(varbio_female_50to60) sum_varbio_female_50to60,
sum(varbio_female_30to50) sum_varbio_female_30to50,
sum(varbio_female_30to60) sum_varbio_female_30to60,
sum(varbio_female_immature) sum_varbio_female_immature
from cb_variancebio_sizegroup
where district = 'East 166'
group by survey_year
order by survey_year;

drop table e166cb_pop_sizegroup_cv;

create table e166cb_pop_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_pop_male_le30 <> 0
	 then ((sqrt(sum_varpop_male_le30))/sum_pop_male_le30)
	 else 0
	 end) cv_pop_male_le30,
(CASE
	 when sum_pop_male_le40 <> 0
	 then ((sqrt(sum_varpop_male_le40))/sum_pop_male_le40)
	 else 0
	 end) cv_pop_male_le40,   
(CASE
	 when sum_pop_male_le50 <> 0
	 then ((sqrt(sum_varpop_male_le50))/sum_pop_male_le50)
	 else 0
	 end) cv_pop_male_le50,
(CASE
	 when sum_pop_male_le60 <> 0
	 then ((sqrt(sum_varpop_male_le60))/sum_pop_male_le60)
	 else 0
	 end) cv_pop_male_le60,
(CASE
	 when sum_pop_male_30to35 <> 0
	 then ((sqrt(sum_varpop_male_30to35))/sum_pop_male_30to35)
	 else 0
	 end) cv_pop_male_30to35,
(CASE
	 when sum_pop_male_35to40 <> 0
	 then ((sqrt(sum_varpop_male_35to40))/sum_pop_male_35to40)
	 else 0
	 end) cv_pop_male_35to40,
(CASE
	 when sum_pop_male_30to40 <> 0
	 then ((sqrt(sum_varpop_male_30to40))/sum_pop_male_30to40)
	 else 0
	 end) cv_pop_male_30to40,
(CASE
	 when sum_pop_male_40to45 <> 0
	 then ((sqrt(sum_varpop_male_40to45))/sum_pop_male_40to45)
	 else 0
	 end) cv_pop_male_40to45,
(CASE
	 when sum_pop_male_45to50 <> 0
	 then ((sqrt(sum_varpop_male_45to50))/sum_pop_male_45to50)
	 else 0
	 end) cv_pop_male_45to50,
(CASE
	 when sum_pop_male_40to50 <> 0
	 then ((sqrt(sum_varpop_male_40to50))/sum_pop_male_40to50)
	 else 0
	 end) cv_pop_male_40to50,
(CASE
	 when sum_pop_male_50to60 <> 0
	 then ((sqrt(sum_varpop_male_50to60))/sum_pop_male_50to60)
	 else 0
	 end) cv_pop_male_50to60,
(CASE
	 when sum_pop_male_30to50 <> 0
	 then ((sqrt(sum_varpop_male_30to50))/sum_pop_male_30to50)
	 else 0
	 end) cv_pop_male_30to50,   
(CASE
	 when sum_pop_male_30to60 <> 0
	 then ((sqrt(sum_varpop_male_30to60))/sum_pop_male_30to60)
	 else 0
	 end) cv_pop_male_30to60,
(CASE
	 when sum_pop_male_le103 <> 0
	 then ((sqrt(sum_varpop_male_le103))/sum_pop_male_le103)
	 else 0
	 end) cv_pop_male_le103,
(CASE
	 when sum_pop_female_le30 <> 0
	 then ((sqrt(sum_varpop_female_le30))/sum_pop_female_le30)
	 else 0
	 end) cv_pop_female_le30,
(CASE
	 when sum_pop_female_le40 <> 0
	 then ((sqrt(sum_varpop_female_le40))/sum_pop_female_le40)
	 else 0
	 end) cv_pop_female_le40,   
(CASE
	 when sum_pop_female_le50 <> 0
	 then ((sqrt(sum_varpop_female_le50))/sum_pop_female_le50)
	 else 0
	 end) cv_pop_female_le50,
(CASE
	 when sum_pop_female_le60 <> 0
	 then ((sqrt(sum_varpop_female_le60))/sum_pop_female_le60)
	 else 0
	 end) cv_pop_female_le60,
(CASE
	 when sum_pop_female_30to35 <> 0
	 then ((sqrt(sum_varpop_female_30to35))/sum_pop_female_30to35)
	 else 0
	 end) cv_pop_female_30to35,
(CASE
	 when sum_pop_female_35to40 <> 0
	 then ((sqrt(sum_varpop_female_35to40))/sum_pop_female_35to40)
	 else 0
	 end) cv_pop_female_35to40,
(CASE
	 when sum_pop_female_30to40 <> 0
	 then ((sqrt(sum_varpop_female_30to40))/sum_pop_female_30to40)
	 else 0
	 end) cv_pop_female_30to40,
(CASE
	 when sum_pop_female_40to45 <> 0
	 then ((sqrt(sum_varpop_female_40to45))/sum_pop_female_40to45)
	 else 0
	 end) cv_pop_female_40to45,
(CASE
	 when sum_pop_female_45to50 <> 0
	 then ((sqrt(sum_varpop_female_45to50))/sum_pop_female_45to50)
	 else 0
	 end) cv_pop_female_45to50,
(CASE
	 when sum_pop_female_40to50 <> 0
	 then ((sqrt(sum_varpop_female_40to50))/sum_pop_female_40to50)
	 else 0
	 end) cv_pop_female_40to50,
(CASE
	 when sum_pop_female_50to60 <> 0
	 then ((sqrt(sum_varpop_female_50to60))/sum_pop_female_50to60)
	 else 0
	 end) cv_pop_female_50to60,
(CASE
	 when sum_pop_female_30to50 <> 0
	 then ((sqrt(sum_varpop_female_30to50))/sum_pop_female_30to50)
	 else 0
	 end) cv_pop_female_30to50,   
(CASE
	 when sum_pop_female_30to60 <> 0
	 then ((sqrt(sum_varpop_female_30to60))/sum_pop_female_30to60)
	 else 0
	 end) cv_pop_female_30to60,
(CASE
	 when sum_pop_female_immature <> 0
	 then ((sqrt(sum_varpop_female_immature))/sum_pop_female_immature)
	 else 0
	 end) cv_pop_female_immature 	 	 	 	  	  	 	 	 	  
from e166cb_varpop_sizegroup_sum a, e166cb_pop_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;	 	 	 

drop table e166cb_bio_sizegroup_cv;

create table e166cb_bio_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_bio_male_le30 <> 0
	 then ((sqrt(sum_varbio_male_le30))/sum_bio_male_le30)
	 else 0
	 end) cv_bio_male_le30,
(CASE
	 when sum_bio_male_le40 <> 0
	 then ((sqrt(sum_varbio_male_le40))/sum_bio_male_le40)
	 else 0
	 end) cv_bio_male_le40,   
(CASE
	 when sum_bio_male_le50 <> 0
	 then ((sqrt(sum_varbio_male_le50))/sum_bio_male_le50)
	 else 0
	 end) cv_bio_male_le50,
(CASE
	 when sum_bio_male_le60 <> 0
	 then ((sqrt(sum_varbio_male_le60))/sum_bio_male_le60)
	 else 0
	 end) cv_bio_male_le60,
(CASE
	 when sum_bio_male_30to35 <> 0
	 then ((sqrt(sum_varbio_male_30to35))/sum_bio_male_30to35)
	 else 0
	 end) cv_bio_male_30to35,
(CASE
	 when sum_bio_male_35to40 <> 0
	 then ((sqrt(sum_varbio_male_35to40))/sum_bio_male_35to40)
	 else 0
	 end) cv_bio_male_35to40,
(CASE
	 when sum_bio_male_30to40 <> 0
	 then ((sqrt(sum_varbio_male_30to40))/sum_bio_male_30to40)
	 else 0
	 end) cv_bio_male_30to40,
(CASE
	 when sum_bio_male_40to45 <> 0
	 then ((sqrt(sum_varbio_male_40to45))/sum_bio_male_40to45)
	 else 0
	 end) cv_bio_male_40to45,
(CASE
	 when sum_bio_male_45to50 <> 0
	 then ((sqrt(sum_varbio_male_45to50))/sum_bio_male_45to50)
	 else 0
	 end) cv_bio_male_45to50,
(CASE
	 when sum_bio_male_40to50 <> 0
	 then ((sqrt(sum_varbio_male_40to50))/sum_bio_male_40to50)
	 else 0
	 end) cv_bio_male_40to50,
(CASE
	 when sum_bio_male_50to60 <> 0
	 then ((sqrt(sum_varbio_male_50to60))/sum_bio_male_50to60)
	 else 0
	 end) cv_bio_male_50to60,
(CASE
	 when sum_bio_male_30to50 <> 0
	 then ((sqrt(sum_varbio_male_30to50))/sum_bio_male_30to50)
	 else 0
	 end) cv_bio_male_30to50,
(CASE
	 when sum_bio_male_30to60 <> 0
	 then ((sqrt(sum_varbio_male_30to60))/sum_bio_male_30to60)
	 else 0
	 end) cv_bio_male_30to60,   
(CASE
	 when sum_bio_male_le103 <> 0
	 then ((sqrt(sum_varbio_male_le103))/sum_bio_male_le103)
	 else 0
	 end) cv_bio_male_le103,
(CASE
	 when sum_bio_female_le30 <> 0
	 then ((sqrt(sum_varbio_female_le30))/sum_bio_female_le30)
	 else 0
	 end) cv_bio_female_le30,
(CASE
	 when sum_bio_female_le40 <> 0
	 then ((sqrt(sum_varbio_female_le40))/sum_bio_female_le40)
	 else 0
	 end) cv_bio_female_le40,   
(CASE
	 when sum_bio_female_le50 <> 0
	 then ((sqrt(sum_varbio_female_le50))/sum_bio_female_le50)
	 else 0
	 end) cv_bio_female_le50,
(CASE
	 when sum_bio_female_le60 <> 0
	 then ((sqrt(sum_varbio_female_le60))/sum_bio_female_le60)
	 else 0
	 end) cv_bio_female_le60,
(CASE
	 when sum_bio_female_30to35 <> 0
	 then ((sqrt(sum_varbio_female_30to35))/sum_bio_female_30to35)
	 else 0
	 end) cv_bio_female_30to35,
(CASE
	 when sum_bio_female_35to40 <> 0
	 then ((sqrt(sum_varbio_female_35to40))/sum_bio_female_35to40)
	 else 0
	 end) cv_bio_female_35to40,
(CASE
	 when sum_bio_female_30to40 <> 0
	 then ((sqrt(sum_varbio_female_30to40))/sum_bio_female_30to40)
	 else 0
	 end) cv_bio_female_30to40,
(CASE
	 when sum_bio_female_40to45 <> 0
	 then ((sqrt(sum_varbio_female_40to45))/sum_bio_female_40to45)
	 else 0
	 end) cv_bio_female_40to45,
(CASE
	 when sum_bio_female_45to50 <> 0
	 then ((sqrt(sum_varbio_female_45to50))/sum_bio_female_45to50)
	 else 0
	 end) cv_bio_female_45to50,
(CASE
	 when sum_bio_female_40to50 <> 0
	 then ((sqrt(sum_varbio_female_40to50))/sum_bio_female_40to50)
	 else 0
	 end) cv_bio_female_40to50,
(CASE
	 when sum_bio_female_50to60 <> 0
	 then ((sqrt(sum_varbio_female_50to60))/sum_bio_female_50to60)
	 else 0
	 end) cv_bio_female_50to60,
(CASE
	 when sum_bio_female_30to50 <> 0
	 then ((sqrt(sum_varbio_female_30to50))/sum_bio_female_30to50)
	 else 0
	 end) cv_bio_female_30to50,
(CASE
	 when sum_bio_female_30to60 <> 0
	 then ((sqrt(sum_varbio_female_30to60))/sum_bio_female_30to60)
	 else 0
	 end) cv_bio_female_30to60,
(CASE
	 when sum_bio_female_immature <> 0
	 then ((sqrt(sum_varbio_female_immature))/sum_bio_female_immature)
	 else 0
	 end) cv_bio_female_immature 
from e166cb_varbio_sizegroup_sum a, e166cb_bio_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;


-- CI calcs

create or replace view e166cb_sizegroup_stderr_pop as
select distinct survey_year,
(sqrt(sum_varpop_male_le30)) stderr_pop_male_le30,
(sqrt(sum_varpop_male_le40)) stderr_pop_male_le40,
(sqrt(sum_varpop_male_le50)) stderr_pop_male_le50,
(sqrt(sum_varpop_male_le60)) stderr_pop_male_le60,
(sqrt(sum_varpop_male_30to35)) stderr_pop_male_30to35,
(sqrt(sum_varpop_male_35to40)) stderr_pop_male_35to40,
(sqrt(sum_varpop_male_30to40)) stderr_pop_male_30to40,
(sqrt(sum_varpop_male_40to45)) stderr_pop_male_40to45,
(sqrt(sum_varpop_male_45to50)) stderr_pop_male_45to50,
(sqrt(sum_varpop_male_40to50)) stderr_pop_male_40to50,
(sqrt(sum_varpop_male_50to60)) stderr_pop_male_50to60,
(sqrt(sum_varpop_male_30to50)) stderr_pop_male_30to50,
(sqrt(sum_varpop_male_30to60)) stderr_pop_male_30to60,
(sqrt(sum_varpop_male_le103)) stderr_pop_male_le103,
(sqrt(sum_varpop_female_le30)) stderr_pop_female_le30,
(sqrt(sum_varpop_female_le40)) stderr_pop_female_le40,
(sqrt(sum_varpop_female_le50)) stderr_pop_female_le50,
(sqrt(sum_varpop_female_le60)) stderr_pop_female_le60,
(sqrt(sum_varpop_female_30to35)) stderr_pop_female_30to35,
(sqrt(sum_varpop_female_35to40)) stderr_pop_female_35to40,
(sqrt(sum_varpop_female_30to40)) stderr_pop_female_30to40,
(sqrt(sum_varpop_female_40to45)) stderr_pop_female_40to45,
(sqrt(sum_varpop_female_45to50)) stderr_pop_female_45to50,
(sqrt(sum_varpop_female_40to50)) stderr_pop_female_40to50,
(sqrt(sum_varpop_female_50to60)) stderr_pop_female_50to60,
(sqrt(sum_varpop_female_30to50)) stderr_pop_female_30to50,
(sqrt(sum_varpop_female_30to60)) stderr_pop_female_30to60,
(sqrt(sum_varpop_female_immature)) stderr_pop_female_immature
from e166cb_varpop_sizegroup_sum;

create or replace view e166cb_sizegroup_stderr_bio as
select distinct survey_year,
(sqrt(sum_varbio_male_le30)) stderr_bio_male_le30,
(sqrt(sum_varbio_male_le40)) stderr_bio_male_le40,
(sqrt(sum_varbio_male_le50)) stderr_bio_male_le50,
(sqrt(sum_varbio_male_le60)) stderr_bio_male_le60,
(sqrt(sum_varbio_male_30to35)) stderr_bio_male_30to35,
(sqrt(sum_varbio_male_35to40)) stderr_bio_male_35to40,
(sqrt(sum_varbio_male_30to40)) stderr_bio_male_30to40,
(sqrt(sum_varbio_male_40to45)) stderr_bio_male_40to45,
(sqrt(sum_varbio_male_45to50)) stderr_bio_male_45to50,
(sqrt(sum_varbio_male_40to50)) stderr_bio_male_40to50,
(sqrt(sum_varbio_male_50to60)) stderr_bio_male_50to60,
(sqrt(sum_varbio_male_30to50)) stderr_bio_male_30to50,
(sqrt(sum_varbio_male_30to60)) stderr_bio_male_30to60,
(sqrt(sum_varbio_male_le103)) stderr_bio_male_le103,
(sqrt(sum_varbio_female_le30)) stderr_bio_female_le30,
(sqrt(sum_varbio_female_le40)) stderr_bio_female_le40,
(sqrt(sum_varbio_female_le50)) stderr_bio_female_le50,
(sqrt(sum_varbio_female_le60)) stderr_bio_female_le60,
(sqrt(sum_varbio_female_30to35)) stderr_bio_female_30to35,
(sqrt(sum_varbio_female_35to40)) stderr_bio_female_35to40,
(sqrt(sum_varbio_female_30to40)) stderr_bio_female_30to40,
(sqrt(sum_varbio_female_40to45)) stderr_bio_female_40to45,
(sqrt(sum_varbio_female_45to50)) stderr_bio_female_45to50,
(sqrt(sum_varbio_female_40to50)) stderr_bio_female_40to50,
(sqrt(sum_varbio_female_50to60)) stderr_bio_female_50to60,
(sqrt(sum_varbio_female_30to50)) stderr_bio_female_30to50,
(sqrt(sum_varbio_female_30to60)) stderr_bio_female_30to60,
(sqrt(sum_varbio_female_immature)) stderr_bio_female_immature
from e166cb_varbio_sizegroup_sum;

drop table e166cb_sizegroup_ci_pop;

create table e166cb_sizegroup_ci_pop as
select distinct survey_year,
(1.96 * stderr_pop_male_le30) ci_pop_male_le30,
(1.96 * stderr_pop_male_le40) ci_pop_male_le40,
(1.96 * stderr_pop_male_le50) ci_pop_male_le50,
(1.96 * stderr_pop_male_le60) ci_pop_male_le60,
(1.96 * stderr_pop_male_30to35) ci_pop_male_30to35,
(1.96 * stderr_pop_male_35to40) ci_pop_male_35to40,
(1.96 * stderr_pop_male_30to40) ci_pop_male_30to40,
(1.96 * stderr_pop_male_40to45) ci_pop_male_40to45,
(1.96 * stderr_pop_male_45to50) ci_pop_male_45to50,
(1.96 * stderr_pop_male_40to50) ci_pop_male_40to50,
(1.96 * stderr_pop_male_50to60) ci_pop_male_50to60,
(1.96 * stderr_pop_male_30to50) ci_pop_male_30to50,
(1.96 * stderr_pop_male_30to60) ci_pop_male_30to60,
(1.96 * stderr_pop_male_le103) ci_pop_male_le103,
(1.96 * stderr_pop_female_le30) ci_pop_female_le30,
(1.96 * stderr_pop_female_le40) ci_pop_female_le40,
(1.96 * stderr_pop_female_le50) ci_pop_female_le50,
(1.96 * stderr_pop_female_le60) ci_pop_female_le60,
(1.96 * stderr_pop_female_30to35) ci_pop_female_30to35,
(1.96 * stderr_pop_female_35to40) ci_pop_female_35to40,
(1.96 * stderr_pop_female_30to40) ci_pop_female_30to40,
(1.96 * stderr_pop_female_40to45) ci_pop_female_40to45,
(1.96 * stderr_pop_female_45to50) ci_pop_female_45to50,
(1.96 * stderr_pop_female_40to50) ci_pop_female_40to50,
(1.96 * stderr_pop_female_50to60) ci_pop_female_50to60,
(1.96 * stderr_pop_female_30to50) ci_pop_female_30to50,
(1.96 * stderr_pop_female_30to60) ci_pop_female_30to60,
(1.96 * stderr_pop_female_immature) ci_pop_female_immature
from e166cb_sizegroup_stderr_pop;

drop table e166cb_sizegroup_ci_bio;

create table e166cb_sizegroup_ci_bio as
select distinct survey_year,
((1.96 * stderr_bio_male_le30)) ci_bio_male_le30,
((1.96 * stderr_bio_male_le40)) ci_bio_male_le40,
((1.96 * stderr_bio_male_le50)) ci_bio_male_le50,
((1.96 * stderr_bio_male_le60)) ci_bio_male_le60,
((1.96 * stderr_bio_male_30to35)) ci_bio_male_30to35,
((1.96 * stderr_bio_male_35to40)) ci_bio_male_35to40,
((1.96 * stderr_bio_male_30to40)) ci_bio_male_30to40,
((1.96 * stderr_bio_male_40to45)) ci_bio_male_40to45,
((1.96 * stderr_bio_male_45to50)) ci_bio_male_45to50,
((1.96 * stderr_bio_male_40to50)) ci_bio_male_40to50,
((1.96 * stderr_bio_male_50to60)) ci_bio_male_50to60,
((1.96 * stderr_bio_male_30to50)) ci_bio_male_30to50,
((1.96 * stderr_bio_male_30to60)) ci_bio_male_30to60,
((1.96 * stderr_bio_male_le103)) ci_bio_male_le103,
((1.96 * stderr_bio_female_le30)) ci_bio_female_le30,
((1.96 * stderr_bio_female_le40)) ci_bio_female_le40,
((1.96 * stderr_bio_female_le50)) ci_bio_female_le50,
((1.96 * stderr_bio_female_le60)) ci_bio_female_le60,
((1.96 * stderr_bio_female_30to35)) ci_bio_female_30to35,
((1.96 * stderr_bio_female_35to40)) ci_bio_female_35to40,
((1.96 * stderr_bio_female_30to40)) ci_bio_female_30to40,
((1.96 * stderr_bio_female_40to45)) ci_bio_female_40to45,
((1.96 * stderr_bio_female_45to50)) ci_bio_female_45to50,
((1.96 * stderr_bio_female_40to50)) ci_bio_female_40to50,
((1.96 * stderr_bio_female_50to60)) ci_bio_female_50to60,
((1.96 * stderr_bio_female_30to50)) ci_bio_female_30to50,
((1.96 * stderr_bio_female_30to60)) ci_bio_female_30to60,
((1.96 * stderr_bio_female_immature)) ci_bio_female_immature
from e166cb_sizegroup_stderr_bio;


--  west of 166W

drop table w166cb_pop_sizegroup;

create table w166cb_pop_sizegroup as
select survey_year,
sum(pop_male_le30) sum_pop_male_le30,
sum(pop_male_le40) sum_pop_male_le40,
sum(pop_male_le50) sum_pop_male_le50,
sum(pop_male_le60) sum_pop_male_le60,
sum(pop_male_30to35) sum_pop_male_30to35,
sum(pop_male_35to40) sum_pop_male_35to40,
sum(pop_male_30to40) sum_pop_male_30to40,
sum(pop_male_40to45) sum_pop_male_40to45,
sum(pop_male_45to50) sum_pop_male_45to50,
sum(pop_male_40to50) sum_pop_male_40to50,
sum(pop_male_50to60) sum_pop_male_50to60,
sum(pop_male_30to50) sum_pop_male_30to50,
sum(pop_male_30to60) sum_pop_male_30to60,
sum(pop_male_le103) sum_pop_male_le103,
sum(pop_female_le30) sum_pop_female_le30,
sum(pop_female_le40) sum_pop_female_le40,
sum(pop_female_le50) sum_pop_female_le50,
sum(pop_female_le60) sum_pop_female_le60,
sum(pop_female_30to35) sum_pop_female_30to35,
sum(pop_female_35to40) sum_pop_female_35to40,
sum(pop_female_30to40) sum_pop_female_30to40,
sum(pop_female_40to45) sum_pop_female_40to45,
sum(pop_female_45to50) sum_pop_female_45to50,
sum(pop_female_40to50) sum_pop_female_40to50,
sum(pop_female_50to60) sum_pop_female_50to60,
sum(pop_female_30to50) sum_pop_female_30to50,
sum(pop_female_30to60) sum_pop_female_30to60,
sum(pop_female_immature) sum_pop_female_immature
from cb_popbystratum_sizegroup
where district in ('Pribilof MTCA','St. Matthew MTCA','West 166')
group by survey_year
order by survey_year;

drop table w166cb_bio_sizegroup;

create table w166cb_bio_sizegroup as
select survey_year,
sum(bio_male_le30) sum_bio_male_le30,
sum(bio_male_le40) sum_bio_male_le40,
sum(bio_male_le50) sum_bio_male_le50,
sum(bio_male_le60) sum_bio_male_le60,
sum(bio_male_30to35) sum_bio_male_30to35,
sum(bio_male_35to40) sum_bio_male_35to40,
sum(bio_male_30to40) sum_bio_male_30to40,
sum(bio_male_40to45) sum_bio_male_40to45,
sum(bio_male_45to50) sum_bio_male_45to50,
sum(bio_male_40to50) sum_bio_male_40to50,
sum(bio_male_50to60) sum_bio_male_50to60,
sum(bio_male_30to50) sum_bio_male_30to50,
sum(bio_male_30to60) sum_bio_male_30to60,
sum(bio_male_le103) sum_bio_male_le103,
sum(bio_female_le30) sum_bio_female_le30,
sum(bio_female_le40) sum_bio_female_le40,
sum(bio_female_le50) sum_bio_female_le50,
sum(bio_female_le60) sum_bio_female_le60,
sum(bio_female_30to35) sum_bio_female_30to35,
sum(bio_female_35to40) sum_bio_female_35to40,
sum(bio_female_30to40) sum_bio_female_30to40,
sum(bio_female_40to45) sum_bio_female_40to45,
sum(bio_female_45to50) sum_bio_female_45to50,
sum(bio_female_40to50) sum_bio_female_40to50,
sum(bio_female_50to60) sum_bio_female_50to60,
sum(bio_female_30to50) sum_bio_female_30to50,
sum(bio_female_30to60) sum_bio_female_30to60,
sum(bio_female_immature) sum_bio_female_immature
from cb_biobystratum_sizegroup
where district in ('Pribilof MTCA','St. Matthew MTCA','West 166')
group by survey_year
order by survey_year;


drop table w166cb_varpop_sizegroup_sum;

create table w166cb_varpop_sizegroup_sum as
select distinct survey_year,
sum(varpop_male_le30) sum_varpop_male_le30,
sum(varpop_male_le40) sum_varpop_male_le40,
sum(varpop_male_le50) sum_varpop_male_le50,
sum(varpop_male_le60) sum_varpop_male_le60,
sum(varpop_male_30to35) sum_varpop_male_30to35,
sum(varpop_male_35to40) sum_varpop_male_35to40,
sum(varpop_male_30to40) sum_varpop_male_30to40,
sum(varpop_male_40to45) sum_varpop_male_40to45,
sum(varpop_male_45to50) sum_varpop_male_45to50,
sum(varpop_male_40to50) sum_varpop_male_40to50,
sum(varpop_male_50to60) sum_varpop_male_50to60,
sum(varpop_male_30to50) sum_varpop_male_30to50,
sum(varpop_male_30to60) sum_varpop_male_30to60,
sum(varpop_male_le103) sum_varpop_male_le103,
sum(varpop_female_le30) sum_varpop_female_le30,
sum(varpop_female_le40) sum_varpop_female_le40,
sum(varpop_female_le50) sum_varpop_female_le50,
sum(varpop_female_le60) sum_varpop_female_le60,
sum(varpop_female_30to35) sum_varpop_female_30to35,
sum(varpop_female_35to40) sum_varpop_female_35to40,
sum(varpop_female_30to40) sum_varpop_female_30to40,
sum(varpop_female_40to45) sum_varpop_female_40to45,
sum(varpop_female_45to50) sum_varpop_female_45to50,
sum(varpop_female_40to50) sum_varpop_female_40to50,
sum(varpop_female_50to60) sum_varpop_female_50to60,
sum(varpop_female_30to50) sum_varpop_female_30to50,
sum(varpop_female_30to60) sum_varpop_female_30to60,
sum(varpop_female_immature) sum_varpop_female_immature
from cb_variancepop_sizegroup
where district in ('Pribilof MTCA','St. Matthew MTCA','West 166')
group by survey_year
order by survey_year;

drop table w166cb_varbio_sizegroup_sum;

create table w166cb_varbio_sizegroup_sum as
select distinct survey_year,
sum(varbio_male_le30) sum_varbio_male_le30,
sum(varbio_male_le40) sum_varbio_male_le40,
sum(varbio_male_le50) sum_varbio_male_le50,
sum(varbio_male_le60) sum_varbio_male_le60,
sum(varbio_male_30to35) sum_varbio_male_30to35,
sum(varbio_male_35to40) sum_varbio_male_35to40,
sum(varbio_male_30to40) sum_varbio_male_30to40,
sum(varbio_male_40to45) sum_varbio_male_40to45,
sum(varbio_male_45to50) sum_varbio_male_45to50,
sum(varbio_male_40to50) sum_varbio_male_40to50,
sum(varbio_male_50to60) sum_varbio_male_50to60,
sum(varbio_male_30to50) sum_varbio_male_30to50,
sum(varbio_male_30to60) sum_varbio_male_30to60,
sum(varbio_male_le103) sum_varbio_male_le103,
sum(varbio_female_le30) sum_varbio_female_le30,
sum(varbio_female_le40) sum_varbio_female_le40,
sum(varbio_female_le50) sum_varbio_female_le50,
sum(varbio_female_le60) sum_varbio_female_le60,
sum(varbio_female_30to35) sum_varbio_female_30to35,
sum(varbio_female_35to40) sum_varbio_female_35to40,
sum(varbio_female_30to40) sum_varbio_female_30to40,
sum(varbio_female_40to45) sum_varbio_female_40to45,
sum(varbio_female_45to50) sum_varbio_female_45to50,
sum(varbio_female_40to50) sum_varbio_female_40to50,
sum(varbio_female_50to60) sum_varbio_female_50to60,
sum(varbio_female_30to50) sum_varbio_female_30to50,
sum(varbio_female_30to60) sum_varbio_female_30to60,
sum(varbio_female_immature) sum_varbio_female_immature
from cb_variancebio_sizegroup
where district in ('Pribilof MTCA','St. Matthew MTCA','West 166')
group by survey_year
order by survey_year;

drop table w166cb_pop_sizegroup_cv;

create table w166cb_pop_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_pop_male_le30 <> 0
	 then ((sqrt(sum_varpop_male_le30))/sum_pop_male_le30)
	 else 0
	 end) cv_pop_male_le30,
(CASE
	 when sum_pop_male_le40 <> 0
	 then ((sqrt(sum_varpop_male_le40))/sum_pop_male_le40)
	 else 0
	 end) cv_pop_male_le40,   
(CASE
	 when sum_pop_male_le50 <> 0
	 then ((sqrt(sum_varpop_male_le50))/sum_pop_male_le50)
	 else 0
	 end) cv_pop_male_le50,
(CASE
	 when sum_pop_male_le60 <> 0
	 then ((sqrt(sum_varpop_male_le60))/sum_pop_male_le60)
	 else 0
	 end) cv_pop_male_le60,
(CASE
	 when sum_pop_male_30to35 <> 0
	 then ((sqrt(sum_varpop_male_30to35))/sum_pop_male_30to35)
	 else 0
	 end) cv_pop_male_30to35,
(CASE
	 when sum_pop_male_35to40 <> 0
	 then ((sqrt(sum_varpop_male_35to40))/sum_pop_male_35to40)
	 else 0
	 end) cv_pop_male_35to40,
(CASE
	 when sum_pop_male_30to40 <> 0
	 then ((sqrt(sum_varpop_male_30to40))/sum_pop_male_30to40)
	 else 0
	 end) cv_pop_male_30to40,
(CASE
	 when sum_pop_male_40to45 <> 0
	 then ((sqrt(sum_varpop_male_40to45))/sum_pop_male_40to45)
	 else 0
	 end) cv_pop_male_40to45,
(CASE
	 when sum_pop_male_45to50 <> 0
	 then ((sqrt(sum_varpop_male_45to50))/sum_pop_male_45to50)
	 else 0
	 end) cv_pop_male_45to50,
(CASE
	 when sum_pop_male_40to50 <> 0
	 then ((sqrt(sum_varpop_male_40to50))/sum_pop_male_40to50)
	 else 0
	 end) cv_pop_male_40to50,
(CASE
	 when sum_pop_male_50to60 <> 0
	 then ((sqrt(sum_varpop_male_50to60))/sum_pop_male_50to60)
	 else 0
	 end) cv_pop_male_50to60,
(CASE
	 when sum_pop_male_30to50 <> 0
	 then ((sqrt(sum_varpop_male_30to50))/sum_pop_male_30to50)
	 else 0
	 end) cv_pop_male_30to50,   
(CASE
	 when sum_pop_male_30to60 <> 0
	 then ((sqrt(sum_varpop_male_30to60))/sum_pop_male_30to60)
	 else 0
	 end) cv_pop_male_30to60,
(CASE
	 when sum_pop_male_le103 <> 0
	 then ((sqrt(sum_varpop_male_le103))/sum_pop_male_le103)
	 else 0
	 end) cv_pop_male_le103,
(CASE
	 when sum_pop_female_le30 <> 0
	 then ((sqrt(sum_varpop_female_le30))/sum_pop_female_le30)
	 else 0
	 end) cv_pop_female_le30,
(CASE
	 when sum_pop_female_le40 <> 0
	 then ((sqrt(sum_varpop_female_le40))/sum_pop_female_le40)
	 else 0
	 end) cv_pop_female_le40,   
(CASE
	 when sum_pop_female_le50 <> 0
	 then ((sqrt(sum_varpop_female_le50))/sum_pop_female_le50)
	 else 0
	 end) cv_pop_female_le50,
(CASE
	 when sum_pop_female_le60 <> 0
	 then ((sqrt(sum_varpop_female_le60))/sum_pop_female_le60)
	 else 0
	 end) cv_pop_female_le60,
(CASE
	 when sum_pop_female_30to35 <> 0
	 then ((sqrt(sum_varpop_female_30to35))/sum_pop_female_30to35)
	 else 0
	 end) cv_pop_female_30to35,
(CASE
	 when sum_pop_female_35to40 <> 0
	 then ((sqrt(sum_varpop_female_35to40))/sum_pop_female_35to40)
	 else 0
	 end) cv_pop_female_35to40,
(CASE
	 when sum_pop_female_30to40 <> 0
	 then ((sqrt(sum_varpop_female_30to40))/sum_pop_female_30to40)
	 else 0
	 end) cv_pop_female_30to40,
(CASE
	 when sum_pop_female_40to45 <> 0
	 then ((sqrt(sum_varpop_female_40to45))/sum_pop_female_40to45)
	 else 0
	 end) cv_pop_female_40to45,
(CASE
	 when sum_pop_female_45to50 <> 0
	 then ((sqrt(sum_varpop_female_45to50))/sum_pop_female_45to50)
	 else 0
	 end) cv_pop_female_45to50,
(CASE
	 when sum_pop_female_40to50 <> 0
	 then ((sqrt(sum_varpop_female_40to50))/sum_pop_female_40to50)
	 else 0
	 end) cv_pop_female_40to50,
(CASE
	 when sum_pop_female_50to60 <> 0
	 then ((sqrt(sum_varpop_female_50to60))/sum_pop_female_50to60)
	 else 0
	 end) cv_pop_female_50to60,
(CASE
	 when sum_pop_female_30to50 <> 0
	 then ((sqrt(sum_varpop_female_30to50))/sum_pop_female_30to50)
	 else 0
	 end) cv_pop_female_30to50,   
(CASE
	 when sum_pop_female_30to60 <> 0
	 then ((sqrt(sum_varpop_female_30to60))/sum_pop_female_30to60)
	 else 0
	 end) cv_pop_female_30to60,
(CASE
	 when sum_pop_female_immature <> 0
	 then ((sqrt(sum_varpop_female_immature))/sum_pop_female_immature)
	 else 0
	 end) cv_pop_female_immature	  	 	 	 	  
from w166cb_varpop_sizegroup_sum a, w166cb_pop_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;	 	 	 

drop table w166cb_bio_sizegroup_cv;

create table w166cb_bio_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_bio_male_le30 <> 0
	 then ((sqrt(sum_varbio_male_le30))/sum_bio_male_le30)
	 else 0
	 end) cv_bio_male_le30,
(CASE
	 when sum_bio_male_le40 <> 0
	 then ((sqrt(sum_varbio_male_le40))/sum_bio_male_le40)
	 else 0
	 end) cv_bio_male_le40,   
(CASE
	 when sum_bio_male_le50 <> 0
	 then ((sqrt(sum_varbio_male_le50))/sum_bio_male_le50)
	 else 0
	 end) cv_bio_male_le50,
(CASE
	 when sum_bio_male_le60 <> 0
	 then ((sqrt(sum_varbio_male_le60))/sum_bio_male_le60)
	 else 0
	 end) cv_bio_male_le60,
(CASE
	 when sum_bio_male_30to35 <> 0
	 then ((sqrt(sum_varbio_male_30to35))/sum_bio_male_30to35)
	 else 0
	 end) cv_bio_male_30to35,
(CASE
	 when sum_bio_male_35to40 <> 0
	 then ((sqrt(sum_varbio_male_35to40))/sum_bio_male_35to40)
	 else 0
	 end) cv_bio_male_35to40,
(CASE
	 when sum_bio_male_30to40 <> 0
	 then ((sqrt(sum_varbio_male_30to40))/sum_bio_male_30to40)
	 else 0
	 end) cv_bio_male_30to40,
(CASE
	 when sum_bio_male_40to45 <> 0
	 then ((sqrt(sum_varbio_male_40to45))/sum_bio_male_40to45)
	 else 0
	 end) cv_bio_male_40to45,
(CASE
	 when sum_bio_male_45to50 <> 0
	 then ((sqrt(sum_varbio_male_45to50))/sum_bio_male_45to50)
	 else 0
	 end) cv_bio_male_45to50,
(CASE
	 when sum_bio_male_40to50 <> 0
	 then ((sqrt(sum_varbio_male_40to50))/sum_bio_male_40to50)
	 else 0
	 end) cv_bio_male_40to50,
(CASE
	 when sum_bio_male_50to60 <> 0
	 then ((sqrt(sum_varbio_male_50to60))/sum_bio_male_50to60)
	 else 0
	 end) cv_bio_male_50to60,
(CASE
	 when sum_bio_male_30to50 <> 0
	 then ((sqrt(sum_varbio_male_30to50))/sum_bio_male_30to50)
	 else 0
	 end) cv_bio_male_30to50,
(CASE
	 when sum_bio_male_30to60 <> 0
	 then ((sqrt(sum_varbio_male_30to60))/sum_bio_male_30to60)
	 else 0
	 end) cv_bio_male_30to60,   
(CASE
	 when sum_bio_male_le103 <> 0
	 then ((sqrt(sum_varbio_male_le103))/sum_bio_male_le103)
	 else 0
	 end) cv_bio_male_le103,
(CASE
	 when sum_bio_female_le30 <> 0
	 then ((sqrt(sum_varbio_female_le30))/sum_bio_female_le30)
	 else 0
	 end) cv_bio_female_le30,
(CASE
	 when sum_bio_female_le40 <> 0
	 then ((sqrt(sum_varbio_female_le40))/sum_bio_female_le40)
	 else 0
	 end) cv_bio_female_le40,   
(CASE
	 when sum_bio_female_le50 <> 0
	 then ((sqrt(sum_varbio_female_le50))/sum_bio_female_le50)
	 else 0
	 end) cv_bio_female_le50,
(CASE
	 when sum_bio_female_le60 <> 0
	 then ((sqrt(sum_varbio_female_le60))/sum_bio_female_le60)
	 else 0
	 end) cv_bio_female_le60,
(CASE
	 when sum_bio_female_30to35 <> 0
	 then ((sqrt(sum_varbio_female_30to35))/sum_bio_female_30to35)
	 else 0
	 end) cv_bio_female_30to35,
(CASE
	 when sum_bio_female_35to40 <> 0
	 then ((sqrt(sum_varbio_female_35to40))/sum_bio_female_35to40)
	 else 0
	 end) cv_bio_female_35to40,
(CASE
	 when sum_bio_female_30to40 <> 0
	 then ((sqrt(sum_varbio_female_30to40))/sum_bio_female_30to40)
	 else 0
	 end) cv_bio_female_30to40,
(CASE
	 when sum_bio_female_40to45 <> 0
	 then ((sqrt(sum_varbio_female_40to45))/sum_bio_female_40to45)
	 else 0
	 end) cv_bio_female_40to45,
(CASE
	 when sum_bio_female_45to50 <> 0
	 then ((sqrt(sum_varbio_female_45to50))/sum_bio_female_45to50)
	 else 0
	 end) cv_bio_female_45to50,
(CASE
	 when sum_bio_female_40to50 <> 0
	 then ((sqrt(sum_varbio_female_40to50))/sum_bio_female_40to50)
	 else 0
	 end) cv_bio_female_40to50,
(CASE
	 when sum_bio_female_50to60 <> 0
	 then ((sqrt(sum_varbio_female_50to60))/sum_bio_female_50to60)
	 else 0
	 end) cv_bio_female_50to60,
(CASE
	 when sum_bio_female_30to50 <> 0
	 then ((sqrt(sum_varbio_female_30to50))/sum_bio_female_30to50)
	 else 0
	 end) cv_bio_female_30to50,
(CASE
	 when sum_bio_female_30to60 <> 0
	 then ((sqrt(sum_varbio_female_30to60))/sum_bio_female_30to60)
	 else 0
	 end) cv_bio_female_30to60,
(CASE
	 when sum_bio_female_immature <> 0
	 then ((sqrt(sum_varbio_female_immature))/sum_bio_female_immature)
	 else 0
	 end) cv_bio_female_immature
from w166cb_varbio_sizegroup_sum a, w166cb_bio_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;


-- CI calcs

create or replace view w166cb_sizegroup_stderr_pop as
select distinct survey_year,
(sqrt(sum_varpop_male_le30)) stderr_pop_male_le30,
(sqrt(sum_varpop_male_le40)) stderr_pop_male_le40,
(sqrt(sum_varpop_male_le50)) stderr_pop_male_le50,
(sqrt(sum_varpop_male_le60)) stderr_pop_male_le60,
(sqrt(sum_varpop_male_30to35)) stderr_pop_male_30to35,
(sqrt(sum_varpop_male_35to40)) stderr_pop_male_35to40,
(sqrt(sum_varpop_male_30to40)) stderr_pop_male_30to40,
(sqrt(sum_varpop_male_40to45)) stderr_pop_male_40to45,
(sqrt(sum_varpop_male_45to50)) stderr_pop_male_45to50,
(sqrt(sum_varpop_male_40to50)) stderr_pop_male_40to50,
(sqrt(sum_varpop_male_50to60)) stderr_pop_male_50to60,
(sqrt(sum_varpop_male_30to50)) stderr_pop_male_30to50,
(sqrt(sum_varpop_male_30to60)) stderr_pop_male_30to60,
(sqrt(sum_varpop_male_le103)) stderr_pop_male_le103,
(sqrt(sum_varpop_female_le30)) stderr_pop_female_le30,
(sqrt(sum_varpop_female_le40)) stderr_pop_female_le40,
(sqrt(sum_varpop_female_le50)) stderr_pop_female_le50,
(sqrt(sum_varpop_female_le60)) stderr_pop_female_le60,
(sqrt(sum_varpop_female_30to35)) stderr_pop_female_30to35,
(sqrt(sum_varpop_female_35to40)) stderr_pop_female_35to40,
(sqrt(sum_varpop_female_30to40)) stderr_pop_female_30to40,
(sqrt(sum_varpop_female_40to45)) stderr_pop_female_40to45,
(sqrt(sum_varpop_female_45to50)) stderr_pop_female_45to50,
(sqrt(sum_varpop_female_40to50)) stderr_pop_female_40to50,
(sqrt(sum_varpop_female_50to60)) stderr_pop_female_50to60,
(sqrt(sum_varpop_female_30to50)) stderr_pop_female_30to50,
(sqrt(sum_varpop_female_30to60)) stderr_pop_female_30to60,
(sqrt(sum_varpop_female_immature)) stderr_pop_female_immature
from w166cb_varpop_sizegroup_sum;

create or replace view w166cb_sizegroup_stderr_bio as
select distinct survey_year,
(sqrt(sum_varbio_male_le30)) stderr_bio_male_le30,
(sqrt(sum_varbio_male_le40)) stderr_bio_male_le40,
(sqrt(sum_varbio_male_le50)) stderr_bio_male_le50,
(sqrt(sum_varbio_male_le60)) stderr_bio_male_le60,
(sqrt(sum_varbio_male_30to35)) stderr_bio_male_30to35,
(sqrt(sum_varbio_male_35to40)) stderr_bio_male_35to40,
(sqrt(sum_varbio_male_30to40)) stderr_bio_male_30to40,
(sqrt(sum_varbio_male_40to45)) stderr_bio_male_40to45,
(sqrt(sum_varbio_male_45to50)) stderr_bio_male_45to50,
(sqrt(sum_varbio_male_40to50)) stderr_bio_male_40to50,
(sqrt(sum_varbio_male_50to60)) stderr_bio_male_50to60,
(sqrt(sum_varbio_male_30to50)) stderr_bio_male_30to50,
(sqrt(sum_varbio_male_30to60)) stderr_bio_male_30to60,
(sqrt(sum_varbio_male_le103)) stderr_bio_male_le103,
(sqrt(sum_varbio_female_le30)) stderr_bio_female_le30,
(sqrt(sum_varbio_female_le40)) stderr_bio_female_le40,
(sqrt(sum_varbio_female_le50)) stderr_bio_female_le50,
(sqrt(sum_varbio_female_le60)) stderr_bio_female_le60,
(sqrt(sum_varbio_female_30to35)) stderr_bio_female_30to35,
(sqrt(sum_varbio_female_35to40)) stderr_bio_female_35to40,
(sqrt(sum_varbio_female_30to40)) stderr_bio_female_30to40,
(sqrt(sum_varbio_female_40to45)) stderr_bio_female_40to45,
(sqrt(sum_varbio_female_45to50)) stderr_bio_female_45to50,
(sqrt(sum_varbio_female_40to50)) stderr_bio_female_40to50,
(sqrt(sum_varbio_female_50to60)) stderr_bio_female_50to60,
(sqrt(sum_varbio_female_30to50)) stderr_bio_female_30to50,
(sqrt(sum_varbio_female_30to60)) stderr_bio_female_30to60,
(sqrt(sum_varbio_female_immature)) stderr_bio_female_immature
from w166cb_varbio_sizegroup_sum;

drop table w166cb_sizegroup_ci_pop;

create table w166cb_sizegroup_ci_pop as
select distinct survey_year,
(1.96 * stderr_pop_male_le30) ci_pop_male_le30,
(1.96 * stderr_pop_male_le40) ci_pop_male_le40,
(1.96 * stderr_pop_male_le50) ci_pop_male_le50,
(1.96 * stderr_pop_male_le60) ci_pop_male_le60,
(1.96 * stderr_pop_male_30to35) ci_pop_male_30to35,
(1.96 * stderr_pop_male_35to40) ci_pop_male_35to40,
(1.96 * stderr_pop_male_30to40) ci_pop_male_30to40,
(1.96 * stderr_pop_male_40to45) ci_pop_male_40to45,
(1.96 * stderr_pop_male_45to50) ci_pop_male_45to50,
(1.96 * stderr_pop_male_40to50) ci_pop_male_40to50,
(1.96 * stderr_pop_male_50to60) ci_pop_male_50to60,
(1.96 * stderr_pop_male_30to50) ci_pop_male_30to50,
(1.96 * stderr_pop_male_30to60) ci_pop_male_30to60,
(1.96 * stderr_pop_male_le103) ci_pop_male_le103,
(1.96 * stderr_pop_female_le30) ci_pop_female_le30,
(1.96 * stderr_pop_female_le40) ci_pop_female_le40,
(1.96 * stderr_pop_female_le50) ci_pop_female_le50,
(1.96 * stderr_pop_female_le60) ci_pop_female_le60,
(1.96 * stderr_pop_female_30to35) ci_pop_female_30to35,
(1.96 * stderr_pop_female_35to40) ci_pop_female_35to40,
(1.96 * stderr_pop_female_30to40) ci_pop_female_30to40,
(1.96 * stderr_pop_female_40to45) ci_pop_female_40to45,
(1.96 * stderr_pop_female_45to50) ci_pop_female_45to50,
(1.96 * stderr_pop_female_40to50) ci_pop_female_40to50,
(1.96 * stderr_pop_female_50to60) ci_pop_female_50to60,
(1.96 * stderr_pop_female_30to50) ci_pop_female_30to50,
(1.96 * stderr_pop_female_30to60) ci_pop_female_30to60,
(1.96 * stderr_pop_female_immature) ci_pop_female_immature
from w166cb_sizegroup_stderr_pop;

drop table w166cb_sizegroup_ci_bio;

create table w166cb_sizegroup_ci_bio as
select distinct survey_year,
((1.96 * stderr_bio_male_le30)) ci_bio_male_le30,
((1.96 * stderr_bio_male_le40)) ci_bio_male_le40,
((1.96 * stderr_bio_male_le50)) ci_bio_male_le50,
((1.96 * stderr_bio_male_le60)) ci_bio_male_le60,
((1.96 * stderr_bio_male_30to35)) ci_bio_male_30to35,
((1.96 * stderr_bio_male_35to40)) ci_bio_male_35to40,
((1.96 * stderr_bio_male_30to40)) ci_bio_male_30to40,
((1.96 * stderr_bio_male_40to45)) ci_bio_male_40to45,
((1.96 * stderr_bio_male_45to50)) ci_bio_male_45to50,
((1.96 * stderr_bio_male_40to50)) ci_bio_male_40to50,
((1.96 * stderr_bio_male_50to60)) ci_bio_male_50to60,
((1.96 * stderr_bio_male_30to50)) ci_bio_male_30to50,
((1.96 * stderr_bio_male_30to60)) ci_bio_male_30to60,
((1.96 * stderr_bio_male_le103)) ci_bio_male_le103,
((1.96 * stderr_bio_female_le30)) ci_bio_female_le30,
((1.96 * stderr_bio_female_le40)) ci_bio_female_le40,
((1.96 * stderr_bio_female_le50)) ci_bio_female_le50,
((1.96 * stderr_bio_female_le60)) ci_bio_female_le60,
((1.96 * stderr_bio_female_30to35)) ci_bio_female_30to35,
((1.96 * stderr_bio_female_35to40)) ci_bio_female_35to40,
((1.96 * stderr_bio_female_30to40)) ci_bio_female_30to40,
((1.96 * stderr_bio_female_40to45)) ci_bio_female_40to45,
((1.96 * stderr_bio_female_45to50)) ci_bio_female_45to50,
((1.96 * stderr_bio_female_40to50)) ci_bio_female_40to50,
((1.96 * stderr_bio_female_50to60)) ci_bio_female_50to60,
((1.96 * stderr_bio_female_30to50)) ci_bio_female_30to50,
((1.96 * stderr_bio_female_30to60)) ci_bio_female_30to60,
((1.96 * stderr_bio_female_immature)) ci_bio_female_immature
from w166cb_sizegroup_stderr_bio;


-- East of 166W

drop table cb_e166_bio_juvenile;

create table cb_e166_bio_juvenile as
select a.survey_year,
sum_bio_male_le30 biomass_male_le30,cv_bio_male_le30 cv_biomass_male_le30,ci_bio_male_le30 ci_biomass_male_le30,
sum_bio_male_le40 biomass_male_le40,cv_bio_male_le40 cv_biomass_male_le40,ci_bio_male_le40 ci_biomass_male_le40,
sum_bio_male_le50 biomass_male_le50,cv_bio_male_le50 cv_biomass_male_le50,ci_bio_male_le50 ci_biomass_male_le50,
sum_bio_male_le60 biomass_male_le60,cv_bio_male_le60 cv_biomass_male_le60,ci_bio_male_le60 ci_biomass_male_le60,
sum_bio_male_30to35 biomass_male_30to35,cv_bio_male_30to35 cv_biomass_male_30to35,ci_bio_male_30to35 ci_biomass_male_30to35,
sum_bio_male_35to40 biomass_male_35to40,cv_bio_male_35to40 cv_biomass_male_35to40,ci_bio_male_35to40 ci_biomass_male_35to40,
sum_bio_male_30to40 biomass_male_30to40,cv_bio_male_30to40 cv_biomass_male_30to40,ci_bio_male_30to40 ci_biomass_male_30to40,
sum_bio_male_40to45 biomass_male_40to45,cv_bio_male_40to45 cv_biomass_male_40to45,ci_bio_male_40to45 ci_biomass_male_40to45,
sum_bio_male_45to50 biomass_male_45to50,cv_bio_male_45to50 cv_biomass_male_45to50,ci_bio_male_45to50 ci_biomass_male_45to50,
sum_bio_male_40to50 biomass_male_40to50,cv_bio_male_40to50 cv_biomass_male_40to50,ci_bio_male_40to50 ci_biomass_male_40to50,
sum_bio_male_50to60 biomass_male_50to60,cv_bio_male_50to60 cv_biomass_male_50to60,ci_bio_male_50to60 ci_biomass_male_50to60,
sum_bio_male_30to50 biomass_male_30to50,cv_bio_male_30to50 cv_biomass_male_30to50,ci_bio_male_30to50 ci_biomass_male_30to50,
sum_bio_male_30to60 biomass_male_30to60,cv_bio_male_30to60 cv_biomass_male_30to60,ci_bio_male_30to60 ci_biomass_male_30to60,
sum_bio_male_le103 biomass_male_le103,cv_bio_male_le103 cv_biomass_male_le103,ci_bio_male_le103 ci_biomass_male_le103,
sum_bio_female_le30 biomass_female_le30,cv_bio_female_le30 cv_biomass_female_le30,ci_bio_female_le30 ci_biomass_female_le30,
sum_bio_female_le40 biomass_female_le40,cv_bio_female_le40 cv_biomass_female_le40,ci_bio_female_le40 ci_biomass_female_le40,
sum_bio_female_le50 biomass_female_le50,cv_bio_female_le50 cv_biomass_female_le50,ci_bio_female_le50 ci_biomass_female_le50,
sum_bio_female_le60 biomass_female_le60,cv_bio_female_le60 cv_biomass_female_le60,ci_bio_female_le60 ci_biomass_female_le60,
sum_bio_female_30to35 biomass_female_30to35,cv_bio_female_30to35 cv_biomass_female_30to35,ci_bio_female_30to35 ci_biomass_female_30to35,
sum_bio_female_35to40 biomass_female_35to40,cv_bio_female_35to40 cv_biomass_female_35to40,ci_bio_female_35to40 ci_biomass_female_35to40,
sum_bio_female_30to40 biomass_female_30to40,cv_bio_female_30to40 cv_biomass_female_30to40,ci_bio_female_30to40 ci_biomass_female_30to40,
sum_bio_female_40to45 biomass_female_40to45,cv_bio_female_40to45 cv_biomass_female_40to45,ci_bio_female_40to45 ci_biomass_female_40to45,
sum_bio_female_45to50 biomass_female_45to50,cv_bio_female_45to50 cv_biomass_female_45to50,ci_bio_female_45to50 ci_biomass_female_45to50,
sum_bio_female_40to50 biomass_female_40to50,cv_bio_female_40to50 cv_biomass_female_40to50,ci_bio_female_40to50 ci_biomass_female_40to50,
sum_bio_female_50to60 biomass_female_50to60,cv_bio_female_50to60 cv_biomass_female_50to60,ci_bio_female_50to60 ci_biomass_female_50to60,
sum_bio_female_30to50 biomass_female_30to50,cv_bio_female_30to50 cv_biomass_female_30to50,ci_bio_female_30to50 ci_biomass_female_30to50,
sum_bio_female_30to60 biomass_female_30to60,cv_bio_female_30to60 cv_biomass_female_30to60,ci_bio_female_30to60 ci_biomass_female_30to60,
sum_bio_female_immature biomass_female_immature,cv_bio_female_immature cv_biomass_female_immature,ci_bio_female_immature ci_biomass_female_immature
from e166cb_bio_sizegroup a, e166cb_sizegroup_ci_bio b,e166cb_bio_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;


drop table cb_e166_pop_juvenile;

create table cb_e166_pop_juvenile as
select a.survey_year,
sum_pop_male_le30 num_male_le30,cv_pop_male_le30 cv_num_male_le30,ci_pop_male_le30 ci_num_male_le30,
sum_pop_male_le40 num_male_le40,cv_pop_male_le40 cv_num_male_le40,ci_pop_male_le40 ci_num_male_le40,
sum_pop_male_le50 num_male_le50,cv_pop_male_le50 cv_num_male_le50,ci_pop_male_le50 ci_num_male_le50,
sum_pop_male_le60 num_male_le60,cv_pop_male_le60 cv_num_male_le60,ci_pop_male_le60 ci_num_male_le60,
sum_pop_male_30to35 num_male_30to35,cv_pop_male_30to35 cv_num_male_30to35,ci_pop_male_30to35 ci_num_male_30to35,
sum_pop_male_35to40 num_male_35to40,cv_pop_male_35to40 cv_num_male_35to40,ci_pop_male_35to40 ci_num_male_35to40,
sum_pop_male_30to40 num_male_30to40,cv_pop_male_30to40 cv_num_male_30to40,ci_pop_male_30to40 ci_num_male_30to40,
sum_pop_male_40to45 num_male_40to45,cv_pop_male_40to45 cv_num_male_40to45,ci_pop_male_40to45 ci_num_male_40to45,
sum_pop_male_45to50 num_male_45to50,cv_pop_male_45to50 cv_num_male_45to50,ci_pop_male_45to50 ci_num_male_45to50,
sum_pop_male_40to50 num_male_40to50,cv_pop_male_40to50 cv_num_male_40to50,ci_pop_male_40to50 ci_num_male_40to50,
sum_pop_male_50to60 num_male_50to60,cv_pop_male_50to60 cv_num_male_50to60,ci_pop_male_50to60 ci_num_male_50to60,
sum_pop_male_30to50 num_male_30to50,cv_pop_male_30to50 cv_num_male_30to50,ci_pop_male_30to50 ci_num_male_30to50,
sum_pop_male_30to60 num_male_30to60,cv_pop_male_30to60 cv_num_male_30to60,ci_pop_male_30to60 ci_num_male_30to60,
sum_pop_male_le103 num_male_le103,cv_pop_male_le103 cv_num_male_le103,ci_pop_male_le103 ci_num_male_le103,
sum_pop_female_le30 num_female_le30,cv_pop_female_le30 cv_num_female_le30,ci_pop_female_le30 ci_num_female_le30,
sum_pop_female_le40 num_female_le40,cv_pop_female_le40 cv_num_female_le40,ci_pop_female_le40 ci_num_female_le40,
sum_pop_female_le50 num_female_le50,cv_pop_female_le50 cv_num_female_le50,ci_pop_female_le50 ci_num_female_le50,
sum_pop_female_le60 num_female_le60,cv_pop_female_le60 cv_num_female_le60,ci_pop_female_le60 ci_num_female_le60,
sum_pop_female_30to35 num_female_30to35,cv_pop_female_30to35 cv_num_female_30to35,ci_pop_female_30to35 ci_num_female_30to35,
sum_pop_female_35to40 num_female_35to40,cv_pop_female_35to40 cv_num_female_35to40,ci_pop_female_35to40 ci_num_female_35to40,
sum_pop_female_30to40 num_female_30to40,cv_pop_female_30to40 cv_num_female_30to40,ci_pop_female_30to40 ci_num_female_30to40,
sum_pop_female_40to45 num_female_40to45,cv_pop_female_40to45 cv_num_female_40to45,ci_pop_female_40to45 ci_num_female_40to45,
sum_pop_female_45to50 num_female_45to50,cv_pop_female_45to50 cv_num_female_45to50,ci_pop_female_45to50 ci_num_female_45to50,
sum_pop_female_40to50 num_female_40to50,cv_pop_female_40to50 cv_num_female_40to50,ci_pop_female_40to50 ci_num_female_40to50,
sum_pop_female_50to60 num_female_50to60,cv_pop_female_50to60 cv_num_female_50to60,ci_pop_female_50to60 ci_num_female_50to60,
sum_pop_female_30to50 num_female_30to50,cv_pop_female_30to50 cv_num_female_30to50,ci_pop_female_30to50 ci_num_female_30to50,
sum_pop_female_30to60 num_female_30to60,cv_pop_female_30to60 cv_num_female_30to60,ci_pop_female_30to60 ci_num_female_30to60,
sum_pop_female_immature num_female_immature,cv_pop_female_immature cv_num_female_immature,ci_pop_female_immature ci_num_female_immature
from e166cb_pop_sizegroup a, e166cb_sizegroup_ci_pop b,e166cb_pop_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;


-- West of 166W

drop table cb_w166_bio_juvenile;

create table cb_w166_bio_juvenile as
select a.survey_year,
sum_bio_male_le30 biomass_male_le30,cv_bio_male_le30 cv_biomass_male_le30,ci_bio_male_le30 ci_biomass_male_le30,
sum_bio_male_le40 biomass_male_le40,cv_bio_male_le40 cv_biomass_male_le40,ci_bio_male_le40 ci_biomass_male_le40,
sum_bio_male_le50 biomass_male_le50,cv_bio_male_le50 cv_biomass_male_le50,ci_bio_male_le50 ci_biomass_male_le50,
sum_bio_male_le60 biomass_male_le60,cv_bio_male_le60 cv_biomass_male_le60,ci_bio_male_le60 ci_biomass_male_le60,
sum_bio_male_30to35 biomass_male_30to35,cv_bio_male_30to35 cv_biomass_male_30to35,ci_bio_male_30to35 ci_biomass_male_30to35,
sum_bio_male_35to40 biomass_male_35to40,cv_bio_male_35to40 cv_biomass_male_35to40,ci_bio_male_35to40 ci_biomass_male_35to40,
sum_bio_male_30to40 biomass_male_30to40,cv_bio_male_30to40 cv_biomass_male_30to40,ci_bio_male_30to40 ci_biomass_male_30to40,
sum_bio_male_40to45 biomass_male_40to45,cv_bio_male_40to45 cv_biomass_male_40to45,ci_bio_male_40to45 ci_biomass_male_40to45,
sum_bio_male_45to50 biomass_male_45to50,cv_bio_male_45to50 cv_biomass_male_45to50,ci_bio_male_45to50 ci_biomass_male_45to50,
sum_bio_male_40to50 biomass_male_40to50,cv_bio_male_40to50 cv_biomass_male_40to50,ci_bio_male_40to50 ci_biomass_male_40to50,
sum_bio_male_50to60 biomass_male_50to60,cv_bio_male_50to60 cv_biomass_male_50to60,ci_bio_male_50to60 ci_biomass_male_50to60,
sum_bio_male_30to50 biomass_male_30to50,cv_bio_male_30to50 cv_biomass_male_30to50,ci_bio_male_30to50 ci_biomass_male_30to50,
sum_bio_male_30to60 biomass_male_30to60,cv_bio_male_30to60 cv_biomass_male_30to60,ci_bio_male_30to60 ci_biomass_male_30to60,
sum_bio_male_le103 biomass_male_le103,cv_bio_male_le103 cv_biomass_male_le103,ci_bio_male_le103 ci_biomass_male_le103,
sum_bio_female_le30 biomass_female_le30,cv_bio_female_le30 cv_biomass_female_le30,ci_bio_female_le30 ci_biomass_female_le30,
sum_bio_female_le40 biomass_female_le40,cv_bio_female_le40 cv_biomass_female_le40,ci_bio_female_le40 ci_biomass_female_le40,
sum_bio_female_le50 biomass_female_le50,cv_bio_female_le50 cv_biomass_female_le50,ci_bio_female_le50 ci_biomass_female_le50,
sum_bio_female_le60 biomass_female_le60,cv_bio_female_le60 cv_biomass_female_le60,ci_bio_female_le60 ci_biomass_female_le60,
sum_bio_female_30to35 biomass_female_30to35,cv_bio_female_30to35 cv_biomass_female_30to35,ci_bio_female_30to35 ci_biomass_female_30to35,
sum_bio_female_35to40 biomass_female_35to40,cv_bio_female_35to40 cv_biomass_female_35to40,ci_bio_female_35to40 ci_biomass_female_35to40,
sum_bio_female_30to40 biomass_female_30to40,cv_bio_female_30to40 cv_biomass_female_30to40,ci_bio_female_30to40 ci_biomass_female_30to40,
sum_bio_female_40to45 biomass_female_40to45,cv_bio_female_40to45 cv_biomass_female_40to45,ci_bio_female_40to45 ci_biomass_female_40to45,
sum_bio_female_45to50 biomass_female_45to50,cv_bio_female_45to50 cv_biomass_female_45to50,ci_bio_female_45to50 ci_biomass_female_45to50,
sum_bio_female_40to50 biomass_female_40to50,cv_bio_female_40to50 cv_biomass_female_40to50,ci_bio_female_40to50 ci_biomass_female_40to50,
sum_bio_female_50to60 biomass_female_50to60,cv_bio_female_50to60 cv_biomass_female_50to60,ci_bio_female_50to60 ci_biomass_female_50to60,
sum_bio_female_30to50 biomass_female_30to50,cv_bio_female_30to50 cv_biomass_female_30to50,ci_bio_female_30to50 ci_biomass_female_30to50,
sum_bio_female_30to60 biomass_female_30to60,cv_bio_female_30to60 cv_biomass_female_30to60,ci_bio_female_30to60 ci_biomass_female_30to60,
sum_bio_female_immature biomass_female_immature,cv_bio_female_immature cv_biomass_female_immature,ci_bio_female_immature ci_biomass_female_immature
from w166cb_bio_sizegroup a, w166cb_sizegroup_ci_bio b,w166cb_bio_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;



drop table cb_w166_pop_juvenile;

create table cb_w166_pop_juvenile as
select a.survey_year,
sum_pop_male_le30 num_male_le30,cv_pop_male_le30 cv_num_male_le30,ci_pop_male_le30 ci_num_male_le30,
sum_pop_male_le40 num_male_le40,cv_pop_male_le40 cv_num_male_le40,ci_pop_male_le40 ci_num_male_le40,
sum_pop_male_le50 num_male_le50,cv_pop_male_le50 cv_num_male_le50,ci_pop_male_le50 ci_num_male_le50,
sum_pop_male_le60 num_male_le60,cv_pop_male_le60 cv_num_male_le60,ci_pop_male_le60 ci_num_male_le60,
sum_pop_male_30to35 num_male_30to35,cv_pop_male_30to35 cv_num_male_30to35,ci_pop_male_30to35 ci_num_male_30to35,
sum_pop_male_35to40 num_male_35to40,cv_pop_male_35to40 cv_num_male_35to40,ci_pop_male_35to40 ci_num_male_35to40,
sum_pop_male_30to40 num_male_30to40,cv_pop_male_30to40 cv_num_male_30to40,ci_pop_male_30to40 ci_num_male_30to40,
sum_pop_male_40to45 num_male_40to45,cv_pop_male_40to45 cv_num_male_40to45,ci_pop_male_40to45 ci_num_male_40to45,
sum_pop_male_45to50 num_male_45to50,cv_pop_male_45to50 cv_num_male_45to50,ci_pop_male_45to50 ci_num_male_45to50,
sum_pop_male_40to50 num_male_40to50,cv_pop_male_40to50 cv_num_male_40to50,ci_pop_male_40to50 ci_num_male_40to50,
sum_pop_male_50to60 num_male_50to60,cv_pop_male_50to60 cv_num_male_50to60,ci_pop_male_50to60 ci_num_male_50to60,
sum_pop_male_30to50 num_male_30to50,cv_pop_male_30to50 cv_num_male_30to50,ci_pop_male_30to50 ci_num_male_30to50,
sum_pop_male_30to60 num_male_30to60,cv_pop_male_30to60 cv_num_male_30to60,ci_pop_male_30to60 ci_num_male_30to60,
sum_pop_male_le103 num_male_le103,cv_pop_male_le103 cv_num_male_le103,ci_pop_male_le103 ci_num_male_le103,
sum_pop_female_le30 num_female_le30,cv_pop_female_le30 cv_num_female_le30,ci_pop_female_le30 ci_num_female_le30,
sum_pop_female_le40 num_female_le40,cv_pop_female_le40 cv_num_female_le40,ci_pop_female_le40 ci_num_female_le40,
sum_pop_female_le50 num_female_le50,cv_pop_female_le50 cv_num_female_le50,ci_pop_female_le50 ci_num_female_le50,
sum_pop_female_le60 num_female_le60,cv_pop_female_le60 cv_num_female_le60,ci_pop_female_le60 ci_num_female_le60,
sum_pop_female_30to35 num_female_30to35,cv_pop_female_30to35 cv_num_female_30to35,ci_pop_female_30to35 ci_num_female_30to35,
sum_pop_female_35to40 num_female_35to40,cv_pop_female_35to40 cv_num_female_35to40,ci_pop_female_35to40 ci_num_female_35to40,
sum_pop_female_30to40 num_female_30to40,cv_pop_female_30to40 cv_num_female_30to40,ci_pop_female_30to40 ci_num_female_30to40,
sum_pop_female_40to45 num_female_40to45,cv_pop_female_40to45 cv_num_female_40to45,ci_pop_female_40to45 ci_num_female_40to45,
sum_pop_female_45to50 num_female_45to50,cv_pop_female_45to50 cv_num_female_45to50,ci_pop_female_45to50 ci_num_female_45to50,
sum_pop_female_40to50 num_female_40to50,cv_pop_female_40to50 cv_num_female_40to50,ci_pop_female_40to50 ci_num_female_40to50,
sum_pop_female_50to60 num_female_50to60,cv_pop_female_50to60 cv_num_female_50to60,ci_pop_female_50to60 ci_num_female_50to60,
sum_pop_female_30to50 num_female_30to50,cv_pop_female_30to50 cv_num_female_30to50,ci_pop_female_30to50 ci_num_female_30to50,
sum_pop_female_30to60 num_female_30to60,cv_pop_female_30to60 cv_num_female_30to60,ci_pop_female_30to60 ci_num_female_30to60,
sum_pop_female_immature num_female_immature,cv_pop_female_immature cv_num_female_immature,ci_pop_female_immature ci_num_female_immature
from w166cb_pop_sizegroup a, w166cb_sizegroup_ci_pop b,w166cb_pop_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;


