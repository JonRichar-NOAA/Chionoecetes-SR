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
/*
drop table cbairdi_mat;

create table cbairdi_mat
(CRAB_ID             NUMBER(12,0),
HAULJOIN             NUMBER(12,0),
CRUISE               NUMBER(6,0),
VESSEL               NUMBER(4,0),
HAUL                 NUMBER(4,0),
STATION              VARCHAR2(10 BYTE),
SPECIES_KFRC         NUMBER(2,0),
SPECIES_CODE         NUMBER(6,0),
SEX                  NUMBER(2,0),
LENGTH               NUMBER(4,1),
WIDTH                NUMBER(4,1),
SHELL_CONDITION      NUMBER(2,0),
EGG_COLOR            NUMBER(2,0),
EGG_CONDITION        NUMBER(2,0),
CLUTCH_SIZE          NUMBER(2,0),
WEIGHT               NUMBER(8,2),
DISEASE_CODE         NUMBER(2,0),
DISEASE_DORSAL       NUMBER(3,0),
DISEASE_VENTRAL      NUMBER(3,0),
DISEASE_LEGS         NUMBER(3,0),
SAMPLING_FACTOR      NUMBER(12,5),
CHELA_HEIGHT         NUMBER(4,1),
MERUS_LENGTH         NUMBER(4,1),
COMMENTS             VARCHAR2(2000 BYTE),
GIS_STATION          VARCHAR2(8 BYTE),
STATION_ID	     VARCHAR2(6 BYTE),
LONGITUDE            FLOAT,
LATITUDE             FLOAT
);


insert into cbairdi_mat
select *
from CBAIRDI_STATION_ALL;
*/


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
		 (trunc(width/1) * 1);


-- Females (done separately from males because need clutch size info)

drop table cb_number_size1_female;

create table cb_number_size1_female as
select c.hauljoin,c.vessel,c.cruise,c.haul,h.gis_station,species_code,
shell_condition, egg_condition, clutch_size,(trunc(width/1) * 1)size1,
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
     egg_condition,
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
select hauljoin,vessel,cruise,haul,gis_station,species_code,shell_condition,
egg_condition, clutch_size,size1,
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
order by cruise,vessel,haul,gis_station,shell_condition, egg_condition, size1;

-- Using actual female maturity in this run, so select for clutch size here
--Note: use of clutch_size>1 excludes BARREN and CLEAN mature females. For same year releases, one would expect
--empty egg cases to be present. For following year releases, by time of survey crab should have extruded clutch
--if it was going too, thus code 1 suggests female crab did not/will not contribute larvae under either circumstance

drop table cb_number_size1_matfem;

create table cb_number_size1_matfem as
select hauljoin, vessel, cruise, haul, gis_station, species_code, size1,
  (sum(CASE
	   			WHEN   clutch_size = 0  
				   THEN number_female_size1
				ELSE 0
				END))  number_female_size1_immature,	

  (sum(CASE
	   			WHEN   clutch_size = 1  
                AND shell_condition <=2
				    THEN number_female_size1
				ELSE 0
				END))  number_female_size1_barren,
  (sum(CASE
	   			WHEN   clutch_size = 1  
                AND shell_condition > 2
				    THEN number_female_size1
				ELSE 0
				END))  number_female_size1_oldnoegg,  
  (sum(CASE
	   			WHEN   clutch_size = 1  
                AND egg_condition = 0
                THEN number_female_size1
				ELSE 0
				END))  number_female_size1_noegg, 
  (sum(CASE
	   			WHEN   clutch_size = 1  
                AND egg_condition = 4
                THEN number_female_size1
				ELSE 0
				END))  number_female_size1_hatched, 
  (sum(CASE
	   			WHEN   clutch_size = 1  
                THEN number_female_size1
				ELSE 0
				END))  number_female_size1_ne_tot,
  (sum(CASE
	   			WHEN   clutch_size > 1  
				   THEN number_female_size1
				ELSE 0
				END))  number_female_size1_ovig,
  (sum(CASE
	   			WHEN   clutch_size = 2  
				   THEN number_female_size1
				ELSE 0
				END))  number_female_size1_trace,
  (sum(CASE
	   			WHEN   clutch_size = 3  
				   THEN number_female_size1
				ELSE 0
				END))  number_female_size1_quarter,
  (sum(CASE
  
	   			WHEN   clutch_size = 4  
				   THEN number_female_size1
				ELSE 0
				END))  number_female_size1_half,
  (sum(CASE
	   			WHEN   clutch_size = 5  
				   THEN number_female_size1
				ELSE 0
				END))  number_female_size1_three_qrt,
  (sum(CASE
	   			WHEN   clutch_size = 6 
				   THEN number_female_size1
				ELSE 0
				END))  number_female_size1_full,
  (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 0
                THEN    number_female_size1
				ELSE 0
				END))  number_female_size1_sc0,
  (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 1
                THEN   number_female_size1
				ELSE 0
				END))  number_female_size1_sc1,
  (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 2
                THEN   number_female_size1
				ELSE 0
				END))  number_female_size1_sc2,
  (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 3
                THEN   number_female_size1
				ELSE 0
				END))  number_female_size1_sc3,
  (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 4
                THEN   number_female_size1
				ELSE 0
				END))  number_female_size1_sc4,
  (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 5
                THEN   number_female_size1
				ELSE 0
				END))  number_female_size1_sc5
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
				END))  wgt_female_size1_immature,	
     (sum(CASE
	   			WHEN   clutch_size = 1  
                AND shell_condition <=2
                THEN wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_barren,
     (sum(CASE
	   			WHEN   clutch_size = 1  
                AND shell_condition > 2
                THEN wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_oldnoegg,  
     (sum(CASE
	   			WHEN   clutch_size = 1  
                AND egg_condition = 0
                THEN wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_noegg, 
     (sum(CASE
	   			WHEN   clutch_size = 1  
                AND egg_condition = 4
                THEN wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_hatched, 
     (sum(CASE
	   			WHEN   clutch_size = 1  
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_ne_tot,
     (sum(CASE
	   			WHEN   clutch_size > 1  
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_ovig,
     (sum(CASE
	   			WHEN   clutch_size = 2  
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_trace,
     (sum(CASE
	   			WHEN   clutch_size = 3  
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_quarter,
     (sum(CASE
	   			WHEN   clutch_size = 4  
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_half,
     (sum(CASE
	   			WHEN   clutch_size = 5  
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_three_qrt,
     (sum(CASE
	   			WHEN   clutch_size = 6  
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_full,
    (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 0
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_sc0,
    (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 1
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_sc1,
    (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 2
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_sc2,
     (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 3
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_sc3,
     (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 4
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_sc4,
     (sum(CASE
	   			WHEN   clutch_size > 0 
                AND    shell_condition = 5
                THEN   wgt_female_size1
				ELSE 0
				END))  wgt_female_size1_sc5
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
  WGT_FEMALE_SIZE1_IMMATURE       NUMBER,
  WGT_FEMALE_SIZE1_BARREN         NUMBER,
  WGT_FEMALE_SIZE1_OLDNOEGG       NUMBER,
  WGT_FEMALE_SIZE1_NOEGG          NUMBER,
  WGT_FEMALE_SIZE1_HATCHED        NUMBER,
  WGT_FEMALE_SIZE1_NE_TOT         NUMBER,
  WGT_FEMALE_SIZE1_OVIG           NUMBER,
  WGT_FEMALE_SIZE1_TRACE          NUMBER,
  WGT_FEMALE_SIZE1_QUARTER        NUMBER,
  WGT_FEMALE_SIZE1_HALF           NUMBER,
  WGT_FEMALE_SIZE1_THREE_QRT      NUMBER,
  WGT_FEMALE_SIZE1_FULL           NUMBER,
  WGT_FEMALE_SIZE1_SC0            NUMBER,
  WGT_FEMALE_SIZE1_SC1            NUMBER,
  WGT_FEMALE_SIZE1_SC2            NUMBER,
  WGT_FEMALE_SIZE1_SC3            NUMBER,
  WGT_FEMALE_SIZE1_SC4            NUMBER,
  WGT_FEMALE_SIZE1_SC5            NUMBER
);
insert into cb_weight_grams_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,shell_condition,
wgt_male_size1,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null
from cb_weight_grams_male;

insert into cb_weight_grams_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,null,
null,wgt_female_size1_immature,wgt_female_size1_barren,wgt_female_size1_oldnoegg,
wgt_female_size1_noegg,wgt_female_size1_hatched,wgt_female_size1_ne_tot,wgt_female_size1_ovig,
wgt_female_size1_trace,wgt_female_size1_quarter,wgt_female_size1_half,wgt_female_size1_three_qrt,
wgt_female_size1_full,wgt_female_size1_sc0,wgt_female_size1_sc1,wgt_female_size1_sc2,wgt_female_size1_sc3,
wgt_female_size1_sc4,wgt_female_size1_sc5
from cb_weight_grams_matfem;

-- convert to metric tons

drop table cb_weight_mt_size1;

create table cb_weight_mt_size1 as
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,shell_condition,
(wgt_male_size1 * 0.000001) mt_male_size1,
(wgt_female_size1_immature * 0.000001) mt_female_size1_immature,
(wgt_female_size1_barren * 0.000001) mt_female_size1_barren,
(wgt_female_size1_oldnoegg * 0.000001) mt_female_size1_oldnoegg,
(wgt_female_size1_noegg * 0.000001) mt_female_size1_noegg,
(wgt_female_size1_hatched * 0.000001) mt_female_size1_hatched,
(wgt_female_size1_ne_tot * 0.000001) mt_female_size1_ne_tot,
(wgt_female_size1_ovig * 0.000001) mt_female_size1_ovig,
(wgt_female_size1_trace * 0.000001) mt_female_size1_trace,
(wgt_female_size1_quarter * 0.000001) mt_female_size1_quarter,
(wgt_female_size1_half * 0.000001) mt_female_size1_half,
(wgt_female_size1_three_qrt * 0.000001) mt_female_size1_three_qrt,
(wgt_female_size1_full * 0.000001) mt_female_size1_full,
(wgt_female_size1_sc0 * 0.000001) mt_female_size1_sc0,
(wgt_female_size1_sc1 * 0.000001) mt_female_size1_sc1,
(wgt_female_size1_sc2 * 0.000001) mt_female_size1_sc2,
(wgt_female_size1_sc3 * 0.000001) mt_female_size1_sc3,
(wgt_female_size1_sc4 * 0.000001) mt_female_size1_sc4,
(wgt_female_size1_sc5 * 0.000001) mt_female_size1_sc5
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
  NUMBER_FEMALE_SIZE1_IMMATURE    NUMBER,
  NUMBER_FEMALE_SIZE1_BARREN      NUMBER,
  NUMBER_FEMALE_SIZE1_OLDNOEGG    NUMBER,
  NUMBER_FEMALE_SIZE1_NOEGG       NUMBER,
  NUMBER_FEMALE_SIZE1_HATCHED     NUMBER,
  NUMBER_FEMALE_SIZE1_NE_TOT      NUMBER,
  NUMBER_FEMALE_SIZE1_OVIG        NUMBER,
  NUMBER_FEMALE_SIZE1_TRACE       NUMBER,
  NUMBER_FEMALE_SIZE1_QUARTER     NUMBER,
  NUMBER_FEMALE_SIZE1_HALF        NUMBER,
  NUMBER_FEMALE_SIZE1_THREE_QRT   NUMBER,
  NUMBER_FEMALE_SIZE1_FULL        NUMBER,
  NUMBER_FEMALE_SIZE1_SC0         NUMBER,
  NUMBER_FEMALE_SIZE1_SC1         NUMBER,
  NUMBER_FEMALE_SIZE1_SC2         NUMBER,
  NUMBER_FEMALE_SIZE1_SC3         NUMBER,
  NUMBER_FEMALE_SIZE1_SC4         NUMBER,
  NUMBER_FEMALE_SIZE1_SC5         NUMBER,
  NUMBER_UNSEXED_SIZE1            NUMBER
);
insert into cb_number_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,shell_condition,
number_male_size1,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null
from cb_number_size1_male;

insert into cb_number_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,null,
null,number_female_size1_immature,number_female_size1_barren,number_female_size1_oldnoegg,
number_female_size1_noegg,number_female_size1_hatched,number_female_size1_ne_tot,
number_female_size1_ovig,number_female_size1_trace,number_female_size1_quarter,number_female_size1_half,
number_female_size1_three_qrt,number_female_size1_full,number_female_size1_sc0,number_female_size1_sc1,
number_female_size1_sc2,number_female_size1_sc3,number_female_size1_sc4,number_female_size1_sc5,null
from cb_number_size1_matfem;

insert into cb_number_size1
select hauljoin,vessel,cruise,haul,gis_station,species_code,size1,null,
null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,number_unsexed_size1
from cb_number_size1_unsexed;

-- This section sums the bairdi Tanner crab catch records by haul, sex,
-- and 1-mm size group.  

drop table cb_number_sizegroup;

create table cb_number_sizegroup as
select hauljoin, vessel, cruise, haul, gis_station, species_code, 
	   
	   (sum(CASE
	   			WHEN  size1 between 0 and 94.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_le94,
     (sum(CASE
	   			WHEN  size1 between 0 and 102.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_le102,	                  
				
	   (sum(CASE
	   			WHEN  size1 between 0 and 109.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_le109,	
				
	   (sum(CASE
	   			WHEN  size1 between 0 and 112.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_le112,	          
	   
	   (sum(CASE
	   			WHEN  size1 between 0 and 119.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_le119,
        
	   (sum(CASE
	   			WHEN  size1 between 103.0 and 124.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_103to124,	        
	   
	   (sum(CASE
	   			WHEN  size1 between 113.0 and 124.9
			    THEN number_male_size1
				ELSE 0
				END))  number_male_113to124,	
		
	   (sum(CASE
	   			WHEN  size1 between 103.0 and 250
			    THEN number_male_size1
				ELSE 0
				END))  number_male_ge103,	                
		
	   (sum(CASE
	   			WHEN  size1 between 113.0 and 250
			    THEN number_male_size1
				ELSE 0
				END))  number_male_ge113,	        
		
	   (sum(CASE
	   			WHEN  size1 between 120.0 and 250
			    THEN number_male_size1
				ELSE 0
				END))  number_male_ge120,	
	   
	   (sum(CASE
	   			WHEN  size1 between 125.0 and 250
			    THEN number_male_size1
				ELSE 0
				END))  number_male_ge125,        
	   
	   (sum(CASE
	   			WHEN  size1 between 110.0 and 250
			    THEN number_male_size1
				ELSE 0
				END))  number_male_ge110,
	   
	   (sum(CASE
	   			WHEN  size1 between 138.0 and 250
			    THEN number_male_size1
				ELSE 0
				END))  number_male_ge138,
	   
	   (sum(CASE
	   			WHEN  size1 between 0 and 250
			    THEN number_male_size1
				ELSE 0
				END))  number_male_total,		
	   sum(number_female_size1_immature) number_female_immature,
     sum(number_female_size1_barren) number_female_barren,
     sum(number_female_size1_oldnoegg) number_female_oldnoegg,
     sum(number_female_size1_noegg) number_female_noegg,
     sum(number_female_size1_hatched) number_female_hatched,
     sum(number_female_size1_ne_tot) number_female_ne_tot,
     sum(number_female_size1_trace) number_female_trace,
     sum(number_female_size1_quarter) number_female_quarter,
     sum(number_female_size1_half) number_female_half,
     sum(number_female_size1_three_qrt) number_female_three_qrt,
     sum(number_female_size1_full) number_female_full,
     sum(number_female_size1_sc0) number_female_sc0,
     sum(number_female_size1_sc1) number_female_sc1,
     sum(number_female_size1_sc2) number_female_sc2,
     sum(number_female_size1_sc3) number_female_sc3,
     sum(number_female_size1_sc4) number_female_sc4,
     sum(number_female_size1_sc5) number_female_sc5,
     sum(number_female_size1_ovig) number_female_ovigerous,
    (sum(number_female_size1_immature)+ sum(number_female_size1_ne_tot)+sum(number_female_size1_ovig)) number_female_total,
	   (sum(CASE
	   			WHEN  size1 between 0 and 250
			    THEN number_unsexed_size1
				ELSE 0
				END))  number_unsexed_total														
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
	   			WHEN  size1 between 0 and 94.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_le94,
	   (sum(CASE
	   			WHEN  size1 between 0 and 102.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_le102,	        
	   (sum(CASE
	   			WHEN  size1 between 0 and 109.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_le109,	
	   (sum(CASE
	   			WHEN  size1 between 0 and 112.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_le112,	          
	   (sum(CASE
	   			WHEN  size1 between 0 and 119.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_le119,
	   (sum(CASE
	   			WHEN  size1 between 103.0 and 124.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_103to124,	        
	   (sum(CASE
	   			WHEN  size1 between 113.0 and 124.9
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_113to124,	
 	   (sum(CASE
	   			WHEN  size1 between 103.0 and 250
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_ge103,               
	   (sum(CASE
	   			WHEN  size1 between 113.0 and 250
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_ge113,        
	   (sum(CASE
	   			WHEN  size1 between 120.0 and 250
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_ge120,
	   (sum(CASE
	   			WHEN  size1 between 125.0 and 250
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_ge125,        
	   (sum(CASE
	   			WHEN  size1 between 110.0 and 250
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_ge110,        
	   (sum(CASE
	   			WHEN  size1 between 138.0 and 250
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_ge138,
	   (sum(CASE
	   			WHEN  size1 between 0 and 250
			    THEN mt_male_size1
				ELSE 0
				END))  mt_male_total,		
	   sum(mt_female_size1_immature) mt_female_immature,
     sum(mt_female_size1_barren) mt_female_barren,
     sum(mt_female_size1_oldnoegg) mt_female_oldnoegg,
     sum(mt_female_size1_noegg) mt_female_noegg,
     sum(mt_female_size1_hatched) mt_female_hatched,
     sum(mt_female_size1_ne_tot) mt_female_ne_tot,
     sum(mt_female_size1_trace) mt_female_trace,
     sum(mt_female_size1_quarter) mt_female_quarter,
     sum(mt_female_size1_half) mt_female_half,
     sum(mt_female_size1_three_qrt) mt_female_three_qrt,
     sum(mt_female_size1_full) mt_female_full,
     sum(mt_female_size1_sc0) mt_female_sc0,
     sum(mt_female_size1_sc1) mt_female_sc1,
     sum(mt_female_size1_sc2) mt_female_sc2,
     sum(mt_female_size1_sc3) mt_female_sc3,
     sum(mt_female_size1_sc4) mt_female_sc4,
     sum(mt_female_size1_sc5) mt_female_sc5,
     sum(mt_female_size1_ovig) mt_female_ovigerous,
     (sum(mt_female_size1_immature)+ sum(mt_female_size1_ne_tot)+ sum(mt_female_size1_ovig)) mt_female_total
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
nvl(number_male_le94,0) number_male_le94,
nvl(number_male_le102,0) number_male_le102,
nvl(number_male_le109,0) number_male_le109,
nvl(number_male_le112,0) number_male_le112,
nvl(number_male_le119,0) number_male_le119,
nvl(number_male_103to124,0) number_male_103to124,
nvl(number_male_113to124,0) number_male_113to124,
nvl(number_male_ge103,0) number_male_ge103,
nvl(number_male_ge113,0) number_male_ge113,
nvl(number_male_ge120,0) number_male_ge120,
nvl(number_male_ge125,0) number_male_ge125,
nvl(number_male_ge110,0) number_male_ge110,
nvl(number_male_ge138,0) number_male_ge138,
nvl(number_male_total,0) number_male_total,
nvl(number_female_immature,0) number_female_immature,
nvl(number_female_barren,0) number_female_barren,
nvl(number_female_oldnoegg,0) number_female_oldnoegg,
nvl(number_female_noegg,0) number_female_noegg,
nvl(number_female_hatched,0) number_female_hatched,
nvl(number_female_ne_tot,0) number_female_ne_tot,
nvl(number_female_trace,0) number_female_trace,
nvl(number_female_quarter,0) number_female_quarter,
nvl(number_female_half,0) number_female_half,
nvl(number_female_three_qrt,0) number_female_three_qrt,
nvl(number_female_full,0) number_female_full,
nvl(number_female_sc0,0) number_female_sc0,
nvl(number_female_sc1,0) number_female_sc1,
nvl(number_female_sc2,0) number_female_sc2,
nvl(number_female_sc3,0) number_female_sc3,
nvl(number_female_sc4,0) number_female_sc4,
nvl(number_female_sc5,0) number_female_sc5,
nvl(number_female_ovigerous,0) number_female_ovigerous,
nvl(number_female_total,0) number_female_total,
nvl(number_unsexed_total,0) number_unsexed_total
from haul_newtimeseries_noretow h full outer join cb_number_sizegroup c
on h.hauljoin = c.hauljoin
where haul_type <> 17;

--  Similarly, by weight.

drop table cb_wgt_sizegroup_union;

create table cb_wgt_sizegroup_union as
select h.hauljoin,h.vessel,h.cruise,h.haul,h.gis_station,survey_year,
nvl(species_code,68560) species_code,
nvl(mt_male_le94,0) mt_male_le94,
nvl(mt_male_le102,0) mt_male_le102,
nvl(mt_male_le109,0) mt_male_le109,
nvl(mt_male_le112,0) mt_male_le112,
nvl(mt_male_le119,0) mt_male_le119,
nvl(mt_male_103to124,0) mt_male_103to124,
nvl(mt_male_113to124,0) mt_male_113to124,
nvl(mt_male_ge103,0) mt_male_ge103,
nvl(mt_male_ge113,0) mt_male_ge113,
nvl(mt_male_ge120,0) mt_male_ge120,
nvl(mt_male_ge125,0) mt_male_ge125,
nvl(mt_male_ge110,0) mt_male_ge110,
nvl(mt_male_ge138,0) mt_male_ge138,
nvl(mt_male_total,0) mt_male_total,
nvl(mt_female_immature,0) mt_female_immature,
nvl(mt_female_barren,0) mt_female_barren,
nvl(mt_female_oldnoegg,0) mt_female_oldnoegg,
nvl(mt_female_noegg,0) mt_female_noegg,
nvl(mt_female_hatched,0) mt_female_hatched,
nvl(mt_female_ne_tot,0) mt_female_ne_tot,
nvl(mt_female_trace,0) mt_female_trace,
nvl(mt_female_quarter,0) mt_female_quarter,
nvl(mt_female_half,0) mt_female_half,
nvl(mt_female_three_qrt,0) mt_female_three_qrt,
nvl(mt_female_full,0) mt_female_full,
nvl(mt_female_sc0,0) mt_female_sc0,
nvl(mt_female_sc1,0) mt_female_sc1,
nvl(mt_female_sc2,0) mt_female_sc2,
nvl(mt_female_sc3,0) mt_female_sc3,
nvl(mt_female_sc4,0) mt_female_sc4,
nvl(mt_female_sc5,0) mt_female_sc5,
nvl(mt_female_ovigerous,0) mt_female_ovigerous,
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
select c.hauljoin,c.vessel,c.cruise,c.haul,c.gis_station,c.survey_year,c.species_code,
(number_male_le94 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_le94,
(number_male_le102 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_le102,
(number_male_le109 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_le109,
(number_male_le112 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_le112,
(number_male_le119 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_le119,
(number_male_103to124 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_103to124,
(number_male_113to124 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_113to124,
(number_male_ge103 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_ge103,
(number_male_ge113 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_ge113,
(number_male_ge120 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_ge120,
(number_male_ge125 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_ge125,
(number_male_ge110 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_ge110,
(number_male_ge138 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_ge138,
(number_male_total / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuenum_total,
(number_female_immature / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_immature,
(number_female_barren / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_barren,
(number_female_oldnoegg / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_oldnoegg,
(number_female_noegg / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_noegg,
(number_female_hatched / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_hatched,
(number_female_ne_tot / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_ne_tot,
(number_female_trace / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_trace,
(number_female_quarter / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_quarter,
(number_female_half / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_half,
(number_female_three_qrt / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_three_qrt,
(number_female_full / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_full,
(number_female_sc0 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_sc0,
(number_female_sc1 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_sc1,
(number_female_sc2 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_sc2,
(number_female_sc3 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_sc3,
(number_female_sc4 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_sc4,
(number_female_sc5 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_sc5,
(number_female_ovigerous / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_ovigerous,
(number_female_total / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuenum_total,
(number_unsexed_total / (((net_width/1000) * distance_fished) * 0.29155335)) unsexed_cpuenum_total
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
(mt_male_le94 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_le94,
(mt_male_le102 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_le102,
(mt_male_le109 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_le109,
(mt_male_le112 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_le112,
(mt_male_le119 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_le119,
(mt_male_103to124 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_103to124,
(mt_male_113to124 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_113to124,
(mt_male_ge103 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_ge103,
(mt_male_ge113 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_ge113,
(mt_male_ge120 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_ge120,
(mt_male_ge125 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_ge125,
(mt_male_ge110 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_ge110,
(mt_male_ge138 / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_ge138,
(mt_male_total / (((net_width/1000) * distance_fished) * 0.29155335)) male_cpuewgt_total,
(mt_female_immature / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_immature,
(mt_female_barren / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_barren,
(mt_female_oldnoegg / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_oldnoegg,
(mt_female_noegg / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_noegg,
(mt_female_hatched / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_hatched,
(mt_female_ne_tot / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_ne_tot,
(mt_female_trace / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_trace,
(mt_female_quarter / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_quarter,
(mt_female_half / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_half,
(mt_female_three_qrt / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_three_qrt,
(mt_female_full / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_full,
(mt_female_sc0 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_sc0,
(mt_female_sc1 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_sc1,
(mt_female_sc2 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_sc2,
(mt_female_sc3 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_sc3,
(mt_female_sc4 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_sc4,
(mt_female_sc5 / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_sc5,
(mt_female_ovigerous / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_ovigerous,
(mt_female_total / (((net_width/1000) * distance_fished) * 0.29155335)) female_cpuewgt_total
from cb_wgt_sizegroup_union c, haul_newtimeseries_noretow h
where c.hauljoin = h.hauljoin
and haul_type <> 17;


drop table cb_meancpuenum_sizegroup;

create table cb_meancpuenum_sizegroup as
select c.survey_year,district,
AVG (male_cpuenum_le94) meancpuenum_male_le94,
AVG (male_cpuenum_le102) meancpuenum_male_le102,
AVG (male_cpuenum_le109) meancpuenum_male_le109,
AVG (male_cpuenum_le112) meancpuenum_male_le112,
AVG (male_cpuenum_le119) meancpuenum_male_le119,
AVG (male_cpuenum_103to124) meancpuenum_male_103to124,
AVG (male_cpuenum_113to124) meancpuenum_male_113to124,
AVG (male_cpuenum_ge103) meancpuenum_male_ge103,
AVG (male_cpuenum_ge113) meancpuenum_male_ge113,
AVG (male_cpuenum_ge120) meancpuenum_male_ge120,
AVG (male_cpuenum_ge125) meancpuenum_male_ge125,
AVG (male_cpuenum_ge110) meancpuenum_male_ge110,
AVG (male_cpuenum_ge138) meancpuenum_male_ge138,
AVG (male_cpuenum_total) meancpuenum_male_total,
AVG (female_cpuenum_immature) meancpuenum_female_immature,
AVG (female_cpuenum_barren) meancpuenum_female_barren,
AVG (female_cpuenum_oldnoegg) meancpuenum_female_oldnoegg,
AVG (female_cpuenum_noegg) meancpuenum_female_noegg,
AVG (female_cpuenum_hatched) meancpuenum_female_hatched,
AVG (female_cpuenum_ne_tot) meancpuenum_female_ne_tot,
AVG (female_cpuenum_trace) meancpuenum_female_trace,
AVG (female_cpuenum_quarter) meancpuenum_female_quarter,
AVG (female_cpuenum_half) meancpuenum_female_half,
AVG (female_cpuenum_three_qrt) meancpuenum_female_three_qrt,
AVG (female_cpuenum_full) meancpuenum_female_full,
AVG (female_cpuenum_sc0) meancpuenum_female_sc0,
AVG (female_cpuenum_sc1) meancpuenum_female_sc1,
AVG (female_cpuenum_sc2) meancpuenum_female_sc2,
AVG (female_cpuenum_sc3) meancpuenum_female_sc3,
AVG (female_cpuenum_sc4) meancpuenum_female_sc4,
AVG (female_cpuenum_sc5) meancpuenum_female_sc5,
AVG (female_cpuenum_ovigerous) meancpuenum_female_ovigerous,
AVG (female_cpuenum_total) meancpuenum_female_total,
AVG (unsexed_cpuenum_total) meancpuenum_unsexed_total,
AVG (male_cpuenum_total + female_cpuenum_total + unsexed_cpuenum_total) meancpuenum_gtotal
from cb_cpuenum_sizegroup c, strata_bairdi_newtimeseries s
where c.gis_station = s.station_id
and c.survey_year = s.survey_year
group by c.survey_year,district;

drop table cb_meancpuewgt_sizegroup;

create table cb_meancpuewgt_sizegroup as
select c.survey_year,district,
AVG (male_cpuewgt_le94) meancpuewgt_male_le94,
AVG (male_cpuewgt_le102) meancpuewgt_male_le102,
AVG (male_cpuewgt_le109) meancpuewgt_male_le109,
AVG (male_cpuewgt_le112) meancpuewgt_male_le112,
AVG (male_cpuewgt_le119) meancpuewgt_male_le119,
AVG (male_cpuewgt_103to124) meancpuewgt_male_103to124,
AVG (male_cpuewgt_113to124) meancpuewgt_male_113to124,
AVG (male_cpuewgt_ge103) meancpuewgt_male_ge103,
AVG (male_cpuewgt_ge113) meancpuewgt_male_ge113,
AVG (male_cpuewgt_ge120) meancpuewgt_male_ge120,
AVG (male_cpuewgt_ge125) meancpuewgt_male_ge125,
AVG (male_cpuewgt_ge110) meancpuewgt_male_ge110,
AVG (male_cpuewgt_ge138) meancpuewgt_male_ge138,
AVG (male_cpuewgt_total) meancpuewgt_male_total,
AVG (female_cpuewgt_immature) meancpuewgt_female_immature,
AVG (female_cpuewgt_barren) meancpuewgt_female_barren,
AVG (female_cpuewgt_oldnoegg) meancpuewgt_female_oldnoegg,
AVG (female_cpuewgt_noegg) meancpuewgt_female_noegg,
AVG (female_cpuewgt_hatched) meancpuewgt_female_hatched,
AVG (female_cpuewgt_ne_tot) meancpuewgt_female_ne_tot,
AVG (female_cpuewgt_trace) meancpuewgt_female_trace,
AVG (female_cpuewgt_quarter) meancpuewgt_female_quarter,
AVG (female_cpuewgt_half) meancpuewgt_female_half,
AVG (female_cpuewgt_three_qrt) meancpuewgt_female_three_qrt,
AVG (female_cpuewgt_full) meancpuewgt_female_full,
AVG (female_cpuewgt_sc0) meancpuewgt_female_sc0,
AVG (female_cpuewgt_sc1) meancpuewgt_female_sc1,
AVG (female_cpuewgt_sc2) meancpuewgt_female_sc2,
AVG (female_cpuewgt_sc3) meancpuewgt_female_sc3,
AVG (female_cpuewgt_sc4) meancpuewgt_female_sc4,
AVG (female_cpuewgt_sc5) meancpuewgt_female_sc5,
AVG (female_cpuewgt_ovigerous) meancpuewgt_female_ovigerous,
AVG (female_cpuewgt_total) meancpuewgt_female_total,
AVG (male_cpuewgt_total + female_cpuewgt_total) meancpuewgt_gtotal
from cb_cpuewgt_sizegroup c, strata_bairdi_newtimeseries s
where c.gis_station = s.station_id
and c.survey_year = s.survey_year
group by c.survey_year,district;


drop table cb_popbystratum_sizegroup;

create table cb_popbystratum_sizegroup as
select distinct c.survey_year,stratum,c.district,
 (meancpuenum_male_le94 * total_area) pop_male_le94,
 (meancpuenum_male_le102 * total_area) pop_male_le102,
 (meancpuenum_male_le109 * total_area) pop_male_le109,
 (meancpuenum_male_le112 * total_area) pop_male_le112,
 (meancpuenum_male_le119 * total_area) pop_male_le119,
 (meancpuenum_male_103to124 * total_area) pop_male_103to124,
 (meancpuenum_male_113to124 * total_area) pop_male_113to124,
 (meancpuenum_male_ge103 * total_area) pop_male_ge103,
 (meancpuenum_male_ge113 * total_area) pop_male_ge113,
 (meancpuenum_male_ge120 * total_area) pop_male_ge120,
 (meancpuenum_male_ge125 * total_area) pop_male_ge125,
 (meancpuenum_male_ge110 * total_area) pop_male_ge110,
 (meancpuenum_male_ge138 * total_area) pop_male_ge138,
 (meancpuenum_male_total * total_area) pop_male_total,
 (meancpuenum_female_immature * total_area) pop_female_immature,
 (meancpuenum_female_barren * total_area) pop_female_barren,
 (meancpuenum_female_oldnoegg * total_area) pop_female_oldnoegg,
 (meancpuenum_female_noegg * total_area) pop_female_noegg,
 (meancpuenum_female_hatched * total_area) pop_female_hatched,
 (meancpuenum_female_ne_tot* total_area) pop_female_ne_tot,
 (meancpuenum_female_trace * total_area) pop_female_trace,
 (meancpuenum_female_quarter * total_area) pop_female_quarter,
 (meancpuenum_female_half * total_area) pop_female_half,
 (meancpuenum_female_three_qrt * total_area) pop_female_three_qrt,
 (meancpuenum_female_full * total_area) pop_female_full,
 (meancpuenum_female_sc0 * total_area) pop_female_sc0,
 (meancpuenum_female_sc1 * total_area) pop_female_sc1,
 (meancpuenum_female_sc2 * total_area) pop_female_sc2,
 (meancpuenum_female_sc3 * total_area) pop_female_sc3,
 (meancpuenum_female_sc4 * total_area) pop_female_sc4,
 (meancpuenum_female_sc5 * total_area) pop_female_sc5,
 (meancpuenum_female_ovigerous * total_area) pop_female_ovigerous,
 (meancpuenum_female_total * total_area) pop_female_total,
 (meancpuenum_unsexed_total * total_area) pop_unsexed_total,
 (meancpuenum_gtotal * total_area) pop_gtotal
from cb_meancpuenum_sizegroup c, strata_bairdi_newtimeseries s
where c.district = s.district
and c.survey_year = s.survey_year
order by survey_year,district;

drop table cb_biobystratum_sizegroup;

create table cb_biobystratum_sizegroup as
select distinct c.survey_year,stratum,c.district,
(meancpuewgt_male_le94 * total_area) bio_male_le94,
(meancpuewgt_male_le102 * total_area) bio_male_le102,
(meancpuewgt_male_le109 * total_area) bio_male_le109,
(meancpuewgt_male_le112 * total_area) bio_male_le112,
(meancpuewgt_male_le119 * total_area) bio_male_le119,
(meancpuewgt_male_103to124 * total_area) bio_male_103to124,
(meancpuewgt_male_113to124 * total_area) bio_male_113to124,
(meancpuewgt_male_ge103 * total_area) bio_male_ge103,
(meancpuewgt_male_ge113 * total_area) bio_male_ge113,
(meancpuewgt_male_ge120 * total_area) bio_male_ge120,
(meancpuewgt_male_ge125 * total_area) bio_male_ge125,
(meancpuewgt_male_ge110 * total_area) bio_male_ge110,
(meancpuewgt_male_ge138 * total_area) bio_male_ge138,
(meancpuewgt_male_total * total_area) bio_male_total,
(meancpuewgt_female_immature * total_area) bio_female_immature,
(meancpuewgt_female_barren * total_area) bio_female_barren,
(meancpuewgt_female_oldnoegg * total_area) bio_female_oldnoegg,
(meancpuewgt_female_noegg * total_area) bio_female_noegg,
(meancpuewgt_female_hatched * total_area) bio_female_hatched,
(meancpuewgt_female_ne_tot * total_area) bio_female_ne_tot,
(meancpuewgt_female_trace * total_area) bio_female_trace,
(meancpuewgt_female_quarter * total_area) bio_female_quarter,
(meancpuewgt_female_half * total_area) bio_female_half,
(meancpuewgt_female_three_qrt * total_area) bio_female_three_qrt,
(meancpuewgt_female_full * total_area) bio_female_full,
(meancpuewgt_female_sc0 * total_area) bio_female_sc0,
(meancpuewgt_female_sc1 * total_area) bio_female_sc1,
(meancpuewgt_female_sc2 * total_area) bio_female_sc2,
(meancpuewgt_female_sc3 * total_area) bio_female_sc3,
(meancpuewgt_female_sc4 * total_area) bio_female_sc4,
(meancpuewgt_female_sc5 * total_area) bio_female_sc5,
(meancpuewgt_female_ovigerous * total_area) bio_female_ovigerous,
(meancpuewgt_female_total * total_area) bio_female_total,
(meancpuewgt_gtotal * total_area) bio_gtotal
from cb_meancpuewgt_sizegroup c, strata_bairdi_newtimeseries s
where c.district = s.district
and c.survey_year = s.survey_year
order by survey_year,district;


drop table cb_varcpuenum_sizegroup;

create table cb_varcpuenum_sizegroup as
select c.survey_year,district,
VARIANCE (male_cpuenum_le94) varcpuenum_male_le94,
VARIANCE (male_cpuenum_le102) varcpuenum_male_le102,
VARIANCE (male_cpuenum_le109) varcpuenum_male_le109,
VARIANCE (male_cpuenum_le112) varcpuenum_male_le112,
VARIANCE (male_cpuenum_le119) varcpuenum_male_le119,
VARIANCE (male_cpuenum_103to124) varcpuenum_male_103to124,
VARIANCE (male_cpuenum_113to124) varcpuenum_male_113to124,
VARIANCE (male_cpuenum_ge103) varcpuenum_male_ge103,
VARIANCE (male_cpuenum_ge113) varcpuenum_male_ge113,
VARIANCE (male_cpuenum_ge120) varcpuenum_male_ge120,
VARIANCE (male_cpuenum_ge125) varcpuenum_male_ge125,
VARIANCE (male_cpuenum_ge110) varcpuenum_male_ge110,
VARIANCE (male_cpuenum_ge138) varcpuenum_male_ge138,
VARIANCE (male_cpuenum_total) varcpuenum_male_total,
VARIANCE (female_cpuenum_immature) varcpuenum_female_immature,
VARIANCE (female_cpuenum_barren) varcpuenum_female_barren,
VARIANCE (female_cpuenum_oldnoegg) varcpuenum_female_oldnoegg,
VARIANCE (female_cpuenum_noegg) varcpuenum_female_noegg,
VARIANCE (female_cpuenum_hatched) varcpuenum_female_hatched,
VARIANCE (female_cpuenum_ne_tot) varcpuenum_female_ne_tot,
VARIANCE (female_cpuenum_trace) varcpuenum_female_trace,
VARIANCE (female_cpuenum_quarter) varcpuenum_female_quarter,
VARIANCE (female_cpuenum_half) varcpuenum_female_half,
VARIANCE (female_cpuenum_three_qrt) varcpuenum_female_three_qrt,
VARIANCE (female_cpuenum_full) varcpuenum_female_full,
VARIANCE (female_cpuenum_sc0) varcpuenum_female_sc0,
VARIANCE (female_cpuenum_sc1) varcpuenum_female_sc1,
VARIANCE (female_cpuenum_sc2) varcpuenum_female_sc2,
VARIANCE (female_cpuenum_sc3) varcpuenum_female_sc3,
VARIANCE (female_cpuenum_sc4) varcpuenum_female_sc4,
VARIANCE (female_cpuenum_sc5) varcpuenum_female_sc5,
VARIANCE (female_cpuenum_ovigerous) varcpuenum_female_ovigerous,
VARIANCE (female_cpuenum_total) varcpuenum_female_total,
VARIANCE (unsexed_cpuenum_total) varcpuenum_unsexed_total,
VARIANCE (male_cpuenum_total + female_cpuenum_total + unsexed_cpuenum_total) varcpuenum_gtotal
from cb_cpuenum_sizegroup c, strata_bairdi_newtimeseries s
where c.gis_station = s.station_id
and c.survey_year = s.survey_year
group by c.survey_year,district;

drop table cb_varcpuewgt_sizegroup;

create table cb_varcpuewgt_sizegroup as
select c.survey_year,district,
VARIANCE (male_cpuewgt_le94) varcpuewgt_male_le94,
VARIANCE (male_cpuewgt_le102) varcpuewgt_male_le102,
VARIANCE (male_cpuewgt_le109) varcpuewgt_male_le109,
VARIANCE (male_cpuewgt_le112) varcpuewgt_male_le112,
VARIANCE (male_cpuewgt_le119) varcpuewgt_male_le119,
VARIANCE (male_cpuewgt_103to124) varcpuewgt_male_103to124,
VARIANCE (male_cpuewgt_113to124) varcpuewgt_male_113to124,
VARIANCE (male_cpuewgt_ge103) varcpuewgt_male_ge103,
VARIANCE (male_cpuewgt_ge113) varcpuewgt_male_ge113,
VARIANCE (male_cpuewgt_ge120) varcpuewgt_male_ge120,
VARIANCE (male_cpuewgt_ge125) varcpuewgt_male_ge125,
VARIANCE (male_cpuewgt_ge110) varcpuewgt_male_ge110,
VARIANCE (male_cpuewgt_ge138) varcpuewgt_male_ge138,
VARIANCE (male_cpuewgt_total) varcpuewgt_male_total,
VARIANCE (female_cpuewgt_immature) varcpuewgt_female_immature,
VARIANCE (female_cpuewgt_barren) varcpuewgt_female_barren,
VARIANCE (female_cpuewgt_oldnoegg) varcpuewgt_female_oldnoegg,
VARIANCE (female_cpuewgt_noegg) varcpuewgt_female_noegg,
VARIANCE (female_cpuewgt_hatched) varcpuewgt_female_hatched,
VARIANCE (female_cpuewgt_ne_tot) varcpuewgt_female_ne_tot,
VARIANCE (female_cpuewgt_trace) varcpuewgt_female_trace,
VARIANCE (female_cpuewgt_quarter) varcpuewgt_female_quarter,
VARIANCE (female_cpuewgt_half) varcpuewgt_female_half,
VARIANCE (female_cpuewgt_three_qrt) varcpuewgt_female_three_qrt,
VARIANCE (female_cpuewgt_full) varcpuewgt_female_full,
VARIANCE (female_cpuewgt_sc0) varcpuewgt_female_sc0,
VARIANCE (female_cpuewgt_sc1) varcpuewgt_female_sc1,
VARIANCE (female_cpuewgt_sc2) varcpuewgt_female_sc2,
VARIANCE (female_cpuewgt_sc3) varcpuewgt_female_sc3,
VARIANCE (female_cpuewgt_sc4) varcpuewgt_female_sc4,
VARIANCE (female_cpuewgt_sc5) varcpuewgt_female_sc5,
VARIANCE (female_cpuewgt_ovigerous) varcpuewgt_female_ovigerous,
VARIANCE (female_cpuewgt_total) varcpuewgt_female_total,
VARIANCE (male_cpuewgt_total + female_cpuewgt_total) varcpuewgt_gtotal
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
((varcpuenum_male_le94 * (power(total_area,2)))/number_tows) varpop_male_le94,
((varcpuenum_male_le102 * (power(total_area,2)))/number_tows) varpop_male_le102,
((varcpuenum_male_le109 * (power(total_area,2)))/number_tows) varpop_male_le109,
((varcpuenum_male_le112 * (power(total_area,2)))/number_tows) varpop_male_le112,
((varcpuenum_male_le119 * (power(total_area,2)))/number_tows) varpop_male_le119,
((varcpuenum_male_103to124 * (power(total_area,2)))/number_tows) varpop_male_103to124,
((varcpuenum_male_113to124 * (power(total_area,2)))/number_tows) varpop_male_113to124,
((varcpuenum_male_ge103 * (power(total_area,2)))/number_tows) varpop_male_ge103,
((varcpuenum_male_ge113 * (power(total_area,2)))/number_tows) varpop_male_ge113,
((varcpuenum_male_ge120 * (power(total_area,2)))/number_tows) varpop_male_ge120,
((varcpuenum_male_ge125 * (power(total_area,2)))/number_tows) varpop_male_ge125,
((varcpuenum_male_ge110 * (power(total_area,2)))/number_tows) varpop_male_ge110,
((varcpuenum_male_ge138 * (power(total_area,2)))/number_tows) varpop_male_ge138,
((varcpuenum_male_total * (power(total_area,2)))/number_tows) varpop_male_total,
((varcpuenum_female_immature * (power(total_area,2)))/number_tows) varpop_female_immature,
((varcpuenum_female_barren * (power(total_area,2)))/number_tows) varpop_female_barren,
((varcpuenum_female_oldnoegg * (power(total_area,2)))/number_tows) varpop_female_oldnoegg,
((varcpuenum_female_noegg * (power(total_area,2)))/number_tows) varpop_female_noegg,
((varcpuenum_female_hatched * (power(total_area,2)))/number_tows) varpop_female_hatched,
((varcpuenum_female_ne_tot * (power(total_area,2)))/number_tows) varpop_female_ne_tot,
((varcpuenum_female_trace * (power(total_area,2)))/number_tows) varpop_female_trace,
((varcpuenum_female_quarter * (power(total_area,2)))/number_tows) varpop_female_quarter,
((varcpuenum_female_half * (power(total_area,2)))/number_tows) varpop_female_half,
((varcpuenum_female_three_qrt * (power(total_area,2)))/number_tows) varpop_female_three_qrt,
((varcpuenum_female_full * (power(total_area,2)))/number_tows) varpop_female_full,
((varcpuenum_female_sc0 * (power(total_area,2)))/number_tows) varpop_female_sc0,
((varcpuenum_female_sc1 * (power(total_area,2)))/number_tows) varpop_female_sc1,
((varcpuenum_female_sc2 * (power(total_area,2)))/number_tows) varpop_female_sc2,
((varcpuenum_female_sc3 * (power(total_area,2)))/number_tows) varpop_female_sc3,
((varcpuenum_female_sc4 * (power(total_area,2)))/number_tows) varpop_female_sc4,
((varcpuenum_female_sc5 * (power(total_area,2)))/number_tows) varpop_female_sc5,
((varcpuenum_female_ovigerous * (power(total_area,2)))/number_tows) varpop_female_ovigerous,
((varcpuenum_female_total * (power(total_area,2)))/number_tows) varpop_female_total,
((varcpuenum_unsexed_total * (power(total_area,2)))/number_tows) varpop_unsexed_total,
((varcpuenum_gtotal * (power(total_area,2)))/number_tows) varpop_gtotal
from strata_bairdi_newtimeseries s, cb_varcpuenum_sizegroup c, cb_haulcount n
where c.district = s.district
and c.district = n.district
and c.survey_year = s.survey_year
and c.survey_year = n.survey_year
order by c.survey_year,stratum;

drop table cb_variancebio_sizegroup;

create table cb_variancebio_sizegroup as
select distinct c.survey_year,stratum,c.district,
((varcpuewgt_male_le94 * (power(total_area,2)))/number_tows) varbio_male_le94,
((varcpuewgt_male_le102 * (power(total_area,2)))/number_tows) varbio_male_le102,
((varcpuewgt_male_le109 * (power(total_area,2)))/number_tows) varbio_male_le109,
((varcpuewgt_male_le112 * (power(total_area,2)))/number_tows) varbio_male_le112,
((varcpuewgt_male_le119 * (power(total_area,2)))/number_tows) varbio_male_le119,
((varcpuewgt_male_103to124 * (power(total_area,2)))/number_tows) varbio_male_103to124,
((varcpuewgt_male_113to124 * (power(total_area,2)))/number_tows) varbio_male_113to124,
((varcpuewgt_male_ge103 * (power(total_area,2)))/number_tows) varbio_male_ge103,
((varcpuewgt_male_ge113 * (power(total_area,2)))/number_tows) varbio_male_ge113,
((varcpuewgt_male_ge120 * (power(total_area,2)))/number_tows) varbio_male_ge120,
((varcpuewgt_male_ge125 * (power(total_area,2)))/number_tows) varbio_male_ge125,
((varcpuewgt_male_ge110 * (power(total_area,2)))/number_tows) varbio_male_ge110,
((varcpuewgt_male_ge138 * (power(total_area,2)))/number_tows) varbio_male_ge138,
((varcpuewgt_male_total * (power(total_area,2)))/number_tows) varbio_male_total,
((varcpuewgt_female_immature * (power(total_area,2)))/number_tows) varbio_female_immature,
((varcpuewgt_female_barren * (power(total_area,2)))/number_tows) varbio_female_barren,
((varcpuewgt_female_oldnoegg * (power(total_area,2)))/number_tows) varbio_female_oldnoegg,
((varcpuewgt_female_noegg * (power(total_area,2)))/number_tows) varbio_female_noegg,
((varcpuewgt_female_hatched * (power(total_area,2)))/number_tows) varbio_female_hatched,
((varcpuewgt_female_ne_tot * (power(total_area,2)))/number_tows) varbio_female_ne_tot,
((varcpuewgt_female_trace * (power(total_area,2)))/number_tows) varbio_female_trace,
((varcpuewgt_female_quarter * (power(total_area,2)))/number_tows) varbio_female_quarter,
((varcpuewgt_female_half * (power(total_area,2)))/number_tows) varbio_female_half,
((varcpuewgt_female_three_qrt * (power(total_area,2)))/number_tows) varbio_female_three_qrt,
((varcpuewgt_female_full * (power(total_area,2)))/number_tows) varbio_female_full,
((varcpuewgt_female_sc0 * (power(total_area,2)))/number_tows) varbio_female_sc0,
((varcpuewgt_female_sc1 * (power(total_area,2)))/number_tows) varbio_female_sc1,
((varcpuewgt_female_sc2 * (power(total_area,2)))/number_tows) varbio_female_sc2,
((varcpuewgt_female_sc3 * (power(total_area,2)))/number_tows) varbio_female_sc3,
((varcpuewgt_female_sc4 * (power(total_area,2)))/number_tows) varbio_female_sc4,
((varcpuewgt_female_sc5 * (power(total_area,2)))/number_tows) varbio_female_sc5,
((varcpuewgt_female_ovigerous * (power(total_area,2)))/number_tows) varbio_female_ovigerous,
((varcpuewgt_female_total * (power(total_area,2)))/number_tows) varbio_female_total,
((varcpuewgt_gtotal * (power(total_area,2)))/number_tows) varbio_gtotal
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
sum(pop_male_le94) sum_pop_male_le94,
sum(pop_male_le102) sum_pop_male_le102,
sum(pop_male_le109) sum_pop_male_le109,
sum(pop_male_le112) sum_pop_male_le112,
sum(pop_male_le119) sum_pop_male_le119,
sum(pop_male_103to124) sum_pop_male_103to124,
sum(pop_male_113to124) sum_pop_male_113to124,
sum(pop_male_ge103) sum_pop_male_ge103,
sum(pop_male_ge113) sum_pop_male_ge113,
sum(pop_male_ge120) sum_pop_male_ge120,
sum(pop_male_ge125) sum_pop_male_ge125,
sum(pop_male_ge110) sum_pop_male_ge110,
sum(pop_male_ge138) sum_pop_male_ge138,
sum(pop_male_total) sum_pop_male_total,
sum(pop_female_immature) sum_pop_female_immature,
sum(pop_female_barren) sum_pop_female_barren,
sum(pop_female_oldnoegg) sum_pop_female_oldnoegg,
sum(pop_female_noegg) sum_pop_female_noegg,
sum(pop_female_hatched) sum_pop_female_hatched,
sum(pop_female_ne_tot) sum_pop_female_ne_tot,
sum(pop_female_trace) sum_pop_female_trace,
sum(pop_female_quarter) sum_pop_female_quarter,
sum(pop_female_half) sum_pop_female_half,
sum(pop_female_three_qrt) sum_pop_female_three_qrt,
sum(pop_female_full) sum_pop_female_full,
sum(pop_female_sc0) sum_pop_female_sc0,
sum(pop_female_sc1) sum_pop_female_sc1,
sum(pop_female_sc2) sum_pop_female_sc2,
sum(pop_female_sc3) sum_pop_female_sc3,
sum(pop_female_sc4) sum_pop_female_sc4,
sum(pop_female_sc5) sum_pop_female_sc5,
sum(pop_female_ovigerous) sum_pop_female_ovigerous,
sum(pop_female_total) sum_pop_female_total,
sum(pop_unsexed_total) sum_pop_unsexed_total,
sum(pop_gtotal) sum_pop_gtotal
from cb_popbystratum_sizegroup
group by survey_year
order by survey_year;

drop table cb_bioall_sizegroup;

create table cb_bioall_sizegroup as
select survey_year,
sum(bio_male_le94) sum_bio_male_le94,
sum(bio_male_le102) sum_bio_male_le102,
sum(bio_male_le109) sum_bio_male_le109,
sum(bio_male_le112) sum_bio_male_le112,
sum(bio_male_le119) sum_bio_male_le119,
sum(bio_male_103to124) sum_bio_male_103to124,
sum(bio_male_113to124) sum_bio_male_113to124,
sum(bio_male_ge103) sum_bio_male_ge103,
sum(bio_male_ge113) sum_bio_male_ge113,
sum(bio_male_ge120) sum_bio_male_ge120,
sum(bio_male_ge125) sum_bio_male_ge125,
sum(bio_male_ge110) sum_bio_male_ge110,
sum(bio_male_ge138) sum_bio_male_ge138,
sum(bio_male_total) sum_bio_male_total,
sum(bio_female_immature) sum_bio_female_immature,
sum(bio_female_barren) sum_bio_female_barren,
sum(bio_female_oldnoegg) sum_bio_female_oldnoegg,
sum(bio_female_noegg) sum_bio_female_noegg,
sum(bio_female_hatched) sum_bio_female_hatched,
sum(bio_female_ne_tot) sum_bio_female_ne_tot,
sum(bio_female_trace) sum_bio_female_trace,
sum(bio_female_quarter) sum_bio_female_quarter,
sum(bio_female_half) sum_bio_female_half,
sum(bio_female_three_qrt) sum_bio_female_three_qrt,
sum(bio_female_full) sum_bio_female_full,
sum(bio_female_sc0) sum_bio_female_sc0,
sum(bio_female_sc1) sum_bio_female_sc1,
sum(bio_female_sc2) sum_bio_female_sc2,
sum(bio_female_sc3) sum_bio_female_sc3,
sum(bio_female_sc4) sum_bio_female_sc4,
sum(bio_female_sc5) sum_bio_female_sc5,
sum(bio_female_ovigerous) sum_bio_female_ovigerous,
sum(bio_female_total) sum_bio_female_total,
sum(bio_gtotal) sum_bio_gtotal
from cb_biobystratum_sizegroup
group by survey_year
order by survey_year;


drop table cb_varpop_sizegroup_sum;

create table cb_varpop_sizegroup_sum as
select distinct survey_year,
sum(varpop_male_le94) sum_varpop_male_le94,
sum(varpop_male_le102) sum_varpop_male_le102,
sum(varpop_male_le109) sum_varpop_male_le109,
sum(varpop_male_le112) sum_varpop_male_le112,
sum(varpop_male_le119) sum_varpop_male_le119,
sum(varpop_male_103to124) sum_varpop_male_103to124,
sum(varpop_male_113to124) sum_varpop_male_113to124,
sum(varpop_male_ge103) sum_varpop_male_ge103,
sum(varpop_male_ge113) sum_varpop_male_ge113,
sum(varpop_male_ge120) sum_varpop_male_ge120,
sum(varpop_male_ge125) sum_varpop_male_ge125,
sum(varpop_male_ge110) sum_varpop_male_ge110,
sum(varpop_male_ge138) sum_varpop_male_ge138,
sum(varpop_male_total) sum_varpop_male_total,
sum(varpop_female_immature) sum_varpop_female_immature,
sum(varpop_female_barren) sum_varpop_female_barren,
sum(varpop_female_oldnoegg) sum_varpop_female_oldnoegg,
sum(varpop_female_noegg) sum_varpop_female_noegg,
sum(varpop_female_hatched) sum_varpop_female_hatched,
sum(varpop_female_ne_tot) sum_varpop_female_ne_tot,
sum(varpop_female_trace) sum_varpop_female_trace,
sum(varpop_female_quarter) sum_varpop_female_quarter,
sum(varpop_female_half) sum_varpop_female_half,
sum(varpop_female_three_qrt) sum_varpop_female_three_qrt,
sum(varpop_female_full) sum_varpop_female_full,
sum(varpop_female_sc0) sum_varpop_female_sc0,
sum(varpop_female_sc1) sum_varpop_female_sc1,
sum(varpop_female_sc2) sum_varpop_female_sc2,
sum(varpop_female_sc3) sum_varpop_female_sc3,
sum(varpop_female_sc4) sum_varpop_female_sc4,
sum(varpop_female_sc5) sum_varpop_female_sc5,
sum(varpop_female_ovigerous) sum_varpop_female_ovigerous,
sum(varpop_female_total) sum_varpop_female_total,
sum(varpop_unsexed_total) sum_varpop_unsexed_total,
sum(varpop_gtotal) sum_varpop_gtotal
from cb_variancepop_sizegroup
group by survey_year
order by survey_year;

drop table cb_varbio_sizegroup_sum;

create table cb_varbio_sizegroup_sum as
select distinct survey_year,
sum(varbio_male_le94) sum_varbio_male_le94,
sum(varbio_male_le102) sum_varbio_male_le102,
sum(varbio_male_le109) sum_varbio_male_le109,
sum(varbio_male_le112) sum_varbio_male_le112,
sum(varbio_male_le119) sum_varbio_male_le119,
sum(varbio_male_103to124) sum_varbio_male_103to124,
sum(varbio_male_113to124) sum_varbio_male_113to124,
sum(varbio_male_ge103) sum_varbio_male_ge103,
sum(varbio_male_ge113) sum_varbio_male_ge113,
sum(varbio_male_ge120) sum_varbio_male_ge120,
sum(varbio_male_ge125) sum_varbio_male_ge125,
sum(varbio_male_ge110) sum_varbio_male_ge110,
sum(varbio_male_ge138) sum_varbio_male_ge138,
sum(varbio_male_total) sum_varbio_male_total,
sum(varbio_female_immature) sum_varbio_female_immature,
sum(varbio_female_barren) sum_varbio_female_barren,
sum(varbio_female_oldnoegg) sum_varbio_female_oldnoegg,
sum(varbio_female_noegg) sum_varbio_female_noegg,
sum(varbio_female_hatched) sum_varbio_female_hatched,
sum(varbio_female_ne_tot) sum_varbio_female_ne_tot,
sum(varbio_female_trace) sum_varbio_female_trace,
sum(varbio_female_quarter) sum_varbio_female_quarter,
sum(varbio_female_half) sum_varbio_female_half,
sum(varbio_female_three_qrt) sum_varbio_female_three_qrt,
sum(varbio_female_full) sum_varbio_female_full,
sum(varbio_female_sc0) sum_varbio_female_sc0,
sum(varbio_female_sc1) sum_varbio_female_sc1,
sum(varbio_female_sc2) sum_varbio_female_sc2,
sum(varbio_female_sc3) sum_varbio_female_sc3,
sum(varbio_female_sc4) sum_varbio_female_sc4,
sum(varbio_female_sc5) sum_varbio_female_sc5,
sum(varbio_female_ovigerous) sum_varbio_female_ovigerous,
sum(varbio_female_total) sum_varbio_female_total,
sum(varbio_gtotal) sum_varbio_gtotal
from cb_variancebio_sizegroup
group by survey_year
order by survey_year;

drop table cb_pop_sizegroup_cv;

create table cb_pop_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_pop_male_le94 <> 0
	 then ((sqrt(sum_varpop_male_le94))/sum_pop_male_le94)
	 else 0
	 end) cv_pop_male_le94,
(CASE
	 when sum_pop_male_le102 <> 0
	 then ((sqrt(sum_varpop_male_le102))/sum_pop_male_le102)
	 else 0
	 end) cv_pop_male_le102,   
(CASE
	 when sum_pop_male_le109 <> 0
	 then ((sqrt(sum_varpop_male_le109))/sum_pop_male_le109)
	 else 0
	 end) cv_pop_male_le109,
(CASE
	 when sum_pop_male_le112 <> 0
	 then ((sqrt(sum_varpop_male_le112))/sum_pop_male_le112)
	 else 0
	 end) cv_pop_male_le112,
(CASE
	 when sum_pop_male_le119 <> 0
	 then ((sqrt(sum_varpop_male_le119))/sum_pop_male_le119)
	 else 0
	 end) cv_pop_male_le119,
(CASE
	 when sum_pop_male_103to124 <> 0
	 then ((sqrt(sum_varpop_male_103to124))/sum_pop_male_103to124)
	 else 0
	 end) cv_pop_male_103to124,   
(CASE
	 when sum_pop_male_113to124 <> 0
	 then ((sqrt(sum_varpop_male_113to124))/sum_pop_male_113to124)
	 else 0
	 end) cv_pop_male_113to124,
(CASE
	 when sum_pop_male_ge103 <> 0
	 then ((sqrt(sum_varpop_male_ge103))/sum_pop_male_ge103)
	 else 0
	 end) cv_pop_male_ge103,	       
(CASE
	 when sum_pop_male_ge113 <> 0
	 then ((sqrt(sum_varpop_male_ge113))/sum_pop_male_ge113)
	 else 0
	 end) cv_pop_male_ge113,	    
(CASE
	 when sum_pop_male_ge120 <> 0
	 then ((sqrt(sum_varpop_male_ge120))/sum_pop_male_ge120)
	 else 0
	 end) cv_pop_male_ge120,	 
(CASE
	 when sum_pop_male_ge125 <> 0
	 then ((sqrt(sum_varpop_male_ge125))/sum_pop_male_ge125)
	 else 0
	 end) cv_pop_male_ge125,	 	    
(CASE
	 when sum_pop_male_ge110 <> 0
	 then ((sqrt(sum_varpop_male_ge110))/sum_pop_male_ge110)
	 else 0
	 end) cv_pop_male_ge110,	 	 
(CASE
	 when sum_pop_male_ge138 <> 0
	 then ((sqrt(sum_varpop_male_ge138))/sum_pop_male_ge138)
	 else 0
	 end) cv_pop_male_ge138,
(CASE
	 when sum_pop_male_total <> 0
	 then ((sqrt(sum_varpop_male_total))/sum_pop_male_total)
	 else 0
	 end) cv_pop_male_total,
(CASE
	 when sum_pop_female_immature <> 0
	 then ((sqrt(sum_varpop_female_immature))/sum_pop_female_immature)
	 else 0
	 end) cv_pop_female_immature,	     
(CASE
	 when sum_pop_female_barren <> 0
	 then ((sqrt(sum_varpop_female_barren))/sum_pop_female_barren)
	 else 0
	 end) cv_pop_female_barren,
(CASE
	 when sum_pop_female_oldnoegg <> 0
	 then ((sqrt(sum_varpop_female_oldnoegg))/sum_pop_female_oldnoegg)
	 else 0
	 end) cv_pop_female_oldnoegg,
(CASE
	 when sum_pop_female_noegg <> 0
	 then ((sqrt(sum_varpop_female_noegg))/sum_pop_female_noegg)
	 else 0
	 end) cv_pop_female_noegg,
(CASE
	 when sum_pop_female_hatched <> 0
	 then ((sqrt(sum_varpop_female_hatched))/sum_pop_female_hatched)
	 else 0
	 end) cv_pop_female_hatched,	  
(CASE
	 when sum_pop_female_ne_tot <> 0
	 then ((sqrt(sum_varpop_female_ne_tot))/sum_pop_female_ne_tot)
	 else 0
	 end) cv_pop_female_ne_tot,	  
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
	 when sum_pop_female_three_qrt <> 0
	 then ((sqrt(sum_varpop_female_three_qrt))/sum_pop_female_three_qrt)
	 else 0
	 end) cv_pop_female_three_qrt,
(CASE
	 when sum_pop_female_full <> 0
	 then ((sqrt(sum_varpop_female_full))/sum_pop_female_full)
	 else 0
	 end) cv_pop_female_full,
(CASE
	 when sum_pop_female_sc0 <> 0
	 then ((sqrt(sum_varpop_female_sc0))/sum_pop_female_sc0)
	 else 0
	 end) cv_pop_female_sc0,
(CASE
	 when sum_pop_female_sc1 <> 0
	 then ((sqrt(sum_varpop_female_sc1))/sum_pop_female_sc1)
	 else 0
	 end) cv_pop_female_sc1,
(CASE
	 when sum_pop_female_sc2 <> 0
	 then ((sqrt(sum_varpop_female_sc2))/sum_pop_female_sc2)
	 else 0
	 end) cv_pop_female_sc2,
(CASE
	 when sum_pop_female_sc3 <> 0
	 then ((sqrt(sum_varpop_female_sc3))/sum_pop_female_sc3)
	 else 0
	 end) cv_pop_female_sc3,
(CASE
	 when sum_pop_female_sc4 <> 0
	 then ((sqrt(sum_varpop_female_sc4))/sum_pop_female_sc4)
	 else 0
	 end) cv_pop_female_sc4,
(CASE
	 when sum_pop_female_sc5 <> 0
	 then ((sqrt(sum_varpop_female_sc5))/sum_pop_female_sc5)
	 else 0
	 end) cv_pop_female_sc5,
(CASE
	 when sum_pop_female_ovigerous <> 0
	 then ((sqrt(sum_varpop_female_ovigerous))/sum_pop_female_ovigerous)
	 else 0
	 end) cv_pop_female_ovigerous,	  
(CASE
	 when sum_pop_female_total <> 0
	 then ((sqrt(sum_varpop_female_total))/sum_pop_female_total)
	 else 0
	 end) cv_pop_female_total,
(CASE
	 when sum_pop_unsexed_total <> 0
	 then ((sqrt(sum_varpop_unsexed_total))/sum_pop_unsexed_total)
	 else 0
	 end) cv_pop_unsexed_total,
(CASE
	 when sum_pop_gtotal <> 0
	 then ((sqrt(sum_varpop_gtotal))/sum_pop_gtotal)
	 else 0
	 end) cv_pop_gtotal	  	 	 	 	  
from cb_varpop_sizegroup_sum a, cb_popall_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;	 	 	 

drop table cb_bio_sizegroup_cv;

create table cb_bio_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_bio_male_le94 <> 0
	 then ((sqrt(sum_varbio_male_le94))/sum_bio_male_le94)
	 else 0
	 end) cv_bio_male_le94,
(CASE
	 when sum_bio_male_le102 <> 0
	 then ((sqrt(sum_varbio_male_le102))/sum_bio_male_le102)
	 else 0
	 end) cv_bio_male_le102,   
(CASE
	 when sum_bio_male_le109 <> 0
	 then ((sqrt(sum_varbio_male_le109))/sum_bio_male_le109)
	 else 0
	 end) cv_bio_male_le109,
(CASE
	 when sum_bio_male_le112 <> 0
	 then ((sqrt(sum_varbio_male_le112))/sum_bio_male_le112)
	 else 0
	 end) cv_bio_male_le112,
(CASE
	 when sum_bio_male_le119 <> 0
	 then ((sqrt(sum_varbio_male_le119))/sum_bio_male_le119)
	 else 0
	 end) cv_bio_male_le119,
(CASE
	 when sum_bio_male_103to124 <> 0
	 then ((sqrt(sum_varbio_male_103to124))/sum_bio_male_103to124)
	 else 0
	 end) cv_bio_male_103to124,
(CASE
	 when sum_bio_male_113to124 <> 0
	 then ((sqrt(sum_varbio_male_113to124))/sum_bio_male_113to124)
	 else 0
	 end) cv_bio_male_113to124,   
(CASE
	 when sum_bio_male_ge103 <> 0
	 then ((sqrt(sum_varbio_male_ge103))/sum_bio_male_ge103)
	 else 0
	 end) cv_bio_male_ge103,	       
(CASE
	 when sum_bio_male_ge113 <> 0
	 then ((sqrt(sum_varbio_male_ge113))/sum_bio_male_ge113)
	 else 0
	 end) cv_bio_male_ge113,	    
(CASE
	 when sum_bio_male_ge120 <> 0
	 then ((sqrt(sum_varbio_male_ge120))/sum_bio_male_ge120)
	 else 0
	 end) cv_bio_male_ge120,	 
(CASE
	 when sum_bio_male_ge125 <> 0
	 then ((sqrt(sum_varbio_male_ge125))/sum_bio_male_ge125)
	 else 0
	 end) cv_bio_male_ge125,	 	    
(CASE
	 when sum_bio_male_ge110 <> 0
	 then ((sqrt(sum_varbio_male_ge110))/sum_bio_male_ge110)
	 else 0
	 end) cv_bio_male_ge110,	 	 
(CASE
	 when sum_bio_male_ge138 <> 0
	 then ((sqrt(sum_varbio_male_ge138))/sum_bio_male_ge138)
	 else 0
	 end) cv_bio_male_ge138,
(CASE
	 when sum_bio_male_total <> 0
	 then ((sqrt(sum_varbio_male_total))/sum_bio_male_total)
	 else 0
	 end) cv_bio_male_total,
(CASE
	 when sum_bio_female_immature <> 0
	 then ((sqrt(sum_varbio_female_immature))/sum_bio_female_immature)
	 else 0
	 end) cv_bio_female_immature,	     
(CASE
	 when sum_bio_female_barren <> 0
	 then ((sqrt(sum_varbio_female_barren))/sum_bio_female_barren)
	 else 0
	 end) cv_bio_female_barren,
(CASE
	 when sum_bio_female_oldnoegg <> 0
	 then ((sqrt(sum_varbio_female_oldnoegg))/sum_bio_female_oldnoegg)
	 else 0
	 end) cv_bio_female_oldnoegg,
(CASE
	 when sum_bio_female_noegg <> 0
	 then ((sqrt(sum_varbio_female_noegg))/sum_bio_female_noegg)
	 else 0
	 end) cv_bio_female_noegg,
(CASE
	 when sum_bio_female_hatched <> 0
	 then ((sqrt(sum_varbio_female_hatched))/sum_bio_female_hatched)
	 else 0
	 end) cv_bio_female_hatched,
(CASE
	 when sum_bio_female_ne_tot <> 0
	 then ((sqrt(sum_varbio_female_ne_tot))/sum_bio_female_ne_tot)
	 else 0
	 end) cv_bio_female_ne_tot,
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
	 when sum_bio_female_three_qrt <> 0
	 then ((sqrt(sum_varbio_female_three_qrt))/sum_bio_female_three_qrt)
	 else 0
	 end) cv_bio_female_three_qrt,
(CASE
	 when sum_bio_female_full <> 0
	 then ((sqrt(sum_varbio_female_full))/sum_bio_female_full)
	 else 0
	 end) cv_bio_female_full,
(CASE
	 when sum_bio_female_sc0 <> 0
	 then ((sqrt(sum_varbio_female_sc0))/sum_bio_female_sc0)
	 else 0
	 end) cv_bio_female_sc0,
(CASE
	 when sum_bio_female_sc1 <> 0
	 then ((sqrt(sum_varbio_female_sc1))/sum_bio_female_sc1)
	 else 0
	 end) cv_bio_female_sc1,
(CASE
	 when sum_bio_female_sc2 <> 0
	 then ((sqrt(sum_varbio_female_sc2))/sum_bio_female_sc2)
	 else 0
	 end) cv_bio_female_sc2,
(CASE
	 when sum_bio_female_sc3 <> 0
	 then ((sqrt(sum_varbio_female_sc3))/sum_bio_female_sc3)
	 else 0
	 end) cv_bio_female_sc3,
(CASE
	 when sum_bio_female_sc4 <> 0
	 then ((sqrt(sum_varbio_female_sc4))/sum_bio_female_sc4)
	 else 0
	 end) cv_bio_female_sc4,
(CASE
	 when sum_bio_female_sc5 <> 0
	 then ((sqrt(sum_varbio_female_sc5))/sum_bio_female_sc5)
	 else 0
	 end) cv_bio_female_sc5,
(CASE
	 when sum_bio_female_ovigerous <> 0
	 then ((sqrt(sum_varbio_female_ovigerous))/sum_bio_female_ovigerous)
	 else 0
	 end) cv_bio_female_ovigerous,	
(CASE
	 when sum_bio_female_total <> 0
	 then ((sqrt(sum_varbio_female_total))/sum_bio_female_total)
	 else 0
	 end) cv_bio_female_total,
(CASE
	 when sum_bio_gtotal <> 0
	 then ((sqrt(sum_varbio_gtotal))/sum_bio_gtotal)
	 else 0
	 end) cv_bio_gtotal	  	 	 	 	  
from cb_varbio_sizegroup_sum a, cb_bioall_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;


-- CI calcs

create or replace view cb_sizegroup_stderr_pop as
select distinct survey_year,
(sqrt(sum_varpop_male_le94)) stderr_pop_male_le94,
(sqrt(sum_varpop_male_le102)) stderr_pop_male_le102,
(sqrt(sum_varpop_male_le109)) stderr_pop_male_le109,
(sqrt(sum_varpop_male_le112)) stderr_pop_male_le112,
(sqrt(sum_varpop_male_le119)) stderr_pop_male_le119,
(sqrt(sum_varpop_male_103to124)) stderr_pop_male_103to124,
(sqrt(sum_varpop_male_113to124)) stderr_pop_male_113to124,
(sqrt(sum_varpop_male_ge103)) stderr_pop_male_ge103,
(sqrt(sum_varpop_male_ge113)) stderr_pop_male_ge113,
(sqrt(sum_varpop_male_ge120)) stderr_pop_male_ge120,
(sqrt(sum_varpop_male_ge125)) stderr_pop_male_ge125,
(sqrt(sum_varpop_male_ge110)) stderr_pop_male_ge110,
(sqrt(sum_varpop_male_ge138)) stderr_pop_male_ge138,
(sqrt(sum_varpop_male_total)) stderr_pop_male_total,
(sqrt(sum_varpop_female_immature)) stderr_pop_female_immature,
(sqrt(sum_varpop_female_barren)) stderr_pop_female_barren,
(sqrt(sum_varpop_female_oldnoegg)) stderr_pop_female_oldnoegg,
(sqrt(sum_varpop_female_noegg)) stderr_pop_female_noegg,
(sqrt(sum_varpop_female_hatched)) stderr_pop_female_hatched,
(sqrt(sum_varpop_female_ne_tot)) stderr_pop_female_ne_tot,
(sqrt(sum_varpop_female_trace)) stderr_pop_female_trace,
(sqrt(sum_varpop_female_quarter)) stderr_pop_female_quarter,
(sqrt(sum_varpop_female_half)) stderr_pop_female_half,
(sqrt(sum_varpop_female_three_qrt)) stderr_pop_female_three_qrt,
(sqrt(sum_varpop_female_full)) stderr_pop_female_full,
(sqrt(sum_varpop_female_sc0)) stderr_pop_female_sc0,
(sqrt(sum_varpop_female_sc1)) stderr_pop_female_sc1,
(sqrt(sum_varpop_female_sc2)) stderr_pop_female_sc2,
(sqrt(sum_varpop_female_sc3)) stderr_pop_female_sc3,
(sqrt(sum_varpop_female_sc4)) stderr_pop_female_sc4,
(sqrt(sum_varpop_female_sc5)) stderr_pop_female_sc5,
(sqrt(sum_varpop_female_ovigerous)) stderr_pop_female_ovigerous,
(sqrt(sum_varpop_female_total)) stderr_pop_female_total,
(sqrt(sum_varpop_unsexed_total)) stderr_pop_unsexed_total,
(sqrt(sum_varpop_gtotal)) stderr_pop_gtotal
from cb_varpop_sizegroup_sum;

create or replace view cb_sizegroup_stderr_bio as
select distinct survey_year,
(sqrt(sum_varbio_male_le94)) stderr_bio_male_le94,
(sqrt(sum_varbio_male_le102)) stderr_bio_male_le102,
(sqrt(sum_varbio_male_le109)) stderr_bio_male_le109,
(sqrt(sum_varbio_male_le112)) stderr_bio_male_le112,
(sqrt(sum_varbio_male_le119)) stderr_bio_male_le119,
(sqrt(sum_varbio_male_103to124)) stderr_bio_male_103to124,
(sqrt(sum_varbio_male_113to124)) stderr_bio_male_113to124,
(sqrt(sum_varbio_male_ge103)) stderr_bio_male_ge103,
(sqrt(sum_varbio_male_ge113)) stderr_bio_male_ge113,
(sqrt(sum_varbio_male_ge120)) stderr_bio_male_ge120,
(sqrt(sum_varbio_male_ge125)) stderr_bio_male_ge125,
(sqrt(sum_varbio_male_ge110)) stderr_bio_male_ge110,
(sqrt(sum_varbio_male_ge138)) stderr_bio_male_ge138,
(sqrt(sum_varbio_male_total)) stderr_bio_male_total,
(sqrt(sum_varbio_female_immature)) stderr_bio_female_immature,
(sqrt(sum_varbio_female_barren)) stderr_bio_female_barren,
(sqrt(sum_varbio_female_oldnoegg)) stderr_bio_female_oldnoegg,
(sqrt(sum_varbio_female_noegg)) stderr_bio_female_noegg,
(sqrt(sum_varbio_female_hatched)) stderr_bio_female_hatched,
(sqrt(sum_varbio_female_ne_tot)) stderr_bio_female_ne_tot,
(sqrt(sum_varbio_female_trace)) stderr_bio_female_trace,
(sqrt(sum_varbio_female_quarter)) stderr_bio_female_quarter,
(sqrt(sum_varbio_female_half)) stderr_bio_female_half,
(sqrt(sum_varbio_female_three_qrt)) stderr_bio_female_three_qrt,
(sqrt(sum_varbio_female_full)) stderr_bio_female_full,
(sqrt(sum_varbio_female_sc0)) stderr_bio_female_sc0,
(sqrt(sum_varbio_female_sc1)) stderr_bio_female_sc1,
(sqrt(sum_varbio_female_sc2)) stderr_bio_female_sc2,
(sqrt(sum_varbio_female_sc3)) stderr_bio_female_sc3,
(sqrt(sum_varbio_female_sc4)) stderr_bio_female_sc4,
(sqrt(sum_varbio_female_sc5)) stderr_bio_female_sc5,
(sqrt(sum_varbio_female_ovigerous)) stderr_bio_female_ovigerous,
(sqrt(sum_varbio_female_total)) stderr_bio_female_total,
(sqrt(sum_varbio_gtotal)) stderr_bio_gtotal
from cb_varbio_sizegroup_sum;

drop table cb_sizegroup_confidence_pop;

create table cb_sizegroup_confidence_pop as
select distinct survey_year,
(1.96 * stderr_pop_male_le94) ci_pop_male_le94,
(1.96 * stderr_pop_male_le102) ci_pop_male_le102,
(1.96 * stderr_pop_male_le109) ci_pop_male_le109,
(1.96 * stderr_pop_male_le112) ci_pop_male_le112,
(1.96 * stderr_pop_male_le119) ci_pop_male_le119,
(1.96 * stderr_pop_male_103to124) ci_pop_male_103to124,
(1.96 * stderr_pop_male_113to124) ci_pop_male_113to124,
(1.96 * stderr_pop_male_ge103) ci_pop_male_ge103,
(1.96 * stderr_pop_male_ge113) ci_pop_male_ge113,
(1.96 * stderr_pop_male_ge120) ci_pop_male_ge120,
(1.96 * stderr_pop_male_ge125) ci_pop_male_ge125,
(1.96 * stderr_pop_male_ge110) ci_pop_male_ge110,
(1.96 * stderr_pop_male_ge138) ci_pop_male_ge138,
(1.96 * stderr_pop_male_total) ci_pop_male_total,
(1.96 * stderr_pop_female_immature) ci_pop_female_immature,
(1.96 * stderr_pop_female_barren) ci_pop_female_barren,
(1.96 * stderr_pop_female_oldnoegg) ci_pop_female_oldnoegg,
(1.96 * stderr_pop_female_noegg) ci_pop_female_noegg,
(1.96 * stderr_pop_female_hatched) ci_pop_female_hatched,
(1.96 * stderr_pop_female_ne_tot) ci_pop_female_ne_tot,
(1.96 * stderr_pop_female_trace) ci_pop_female_trace,
(1.96 * stderr_pop_female_quarter) ci_pop_female_quarter,
(1.96 * stderr_pop_female_half) ci_pop_female_half,
(1.96 * stderr_pop_female_three_qrt) ci_pop_female_three_qrt,
(1.96 * stderr_pop_female_full) ci_pop_female_full,
(1.96 * stderr_pop_female_sc0) ci_pop_female_sc0,
(1.96 * stderr_pop_female_sc1) ci_pop_female_sc1,
(1.96 * stderr_pop_female_sc2) ci_pop_female_sc2,
(1.96 * stderr_pop_female_sc3) ci_pop_female_sc3,
(1.96 * stderr_pop_female_sc4) ci_pop_female_sc4,
(1.96 * stderr_pop_female_sc5) ci_pop_female_sc5,
(1.96 * stderr_pop_female_ovigerous) ci_pop_female_ovigerous,
(1.96 * stderr_pop_female_total) ci_pop_female_total,
(1.96 * stderr_pop_unsexed_total) ci_pop_unsexed_total,
(1.96 * stderr_pop_gtotal) ci_pop_gtotal
from cb_sizegroup_stderr_pop;

drop table cb_sizegroup_confidence_bio;

create table cb_sizegroup_confidence_bio as
select distinct survey_year,
((1.96 * stderr_bio_male_le94)) ci_bio_male_le94,
((1.96 * stderr_bio_male_le102)) ci_bio_male_le102,
((1.96 * stderr_bio_male_le109)) ci_bio_male_le109,
((1.96 * stderr_bio_male_le112)) ci_bio_male_le112,
((1.96 * stderr_bio_male_le119)) ci_bio_male_le119,
((1.96 * stderr_bio_male_103to124)) ci_bio_male_103to124,
((1.96 * stderr_bio_male_113to124)) ci_bio_male_113to124,
((1.96 * stderr_bio_male_ge103)) ci_bio_male_ge103,
((1.96 * stderr_bio_male_ge113)) ci_bio_male_ge113,
((1.96 * stderr_bio_male_ge120)) ci_bio_male_ge120,
((1.96 * stderr_bio_male_ge125)) ci_bio_male_ge125,
((1.96 * stderr_bio_male_ge110)) ci_bio_male_ge110,
((1.96 * stderr_bio_male_ge138)) ci_bio_male_ge138,
((1.96 * stderr_bio_male_total)) ci_bio_male_total,
((1.96 * stderr_bio_female_immature)) ci_bio_female_immature,
((1.96 * stderr_bio_female_barren)) ci_bio_female_barren,
((1.96 * stderr_bio_female_oldnoegg)) ci_bio_female_oldnoegg,
((1.96 * stderr_bio_female_noegg)) ci_bio_female_noegg,
((1.96 * stderr_bio_female_hatched)) ci_bio_female_hatched,
((1.96 * stderr_bio_female_ne_tot)) ci_bio_female_ne_tot,
((1.96 * stderr_bio_female_trace)) ci_bio_female_trace,
((1.96 * stderr_bio_female_quarter)) ci_bio_female_quarter,
((1.96 * stderr_bio_female_half)) ci_bio_female_half,
((1.96 * stderr_bio_female_three_qrt)) ci_bio_female_three_qrt,
((1.96 * stderr_bio_female_full)) ci_bio_female_full,
((1.96 * stderr_bio_female_sc0)) ci_bio_female_sc0,
((1.96 * stderr_bio_female_sc1)) ci_bio_female_sc1,
((1.96 * stderr_bio_female_sc2)) ci_bio_female_sc2,
((1.96 * stderr_bio_female_sc3)) ci_bio_female_sc3,
((1.96 * stderr_bio_female_sc4)) ci_bio_female_sc4,
((1.96 * stderr_bio_female_sc5)) ci_bio_female_sc5,
((1.96 * stderr_bio_female_ovigerous)) ci_bio_female_ovigerous,
((1.96 * stderr_bio_female_total)) ci_bio_female_total,
((1.96 * stderr_bio_gtotal)) ci_bio_gtotal
from cb_sizegroup_stderr_bio;


-- Final output for stocks

-- All Districts Combined

drop table cb_ebs_bio_sc_cs_v2024;

create table cb_ebs_bio_sc_cs_v2024 as
select a.survey_year,
sum_bio_male_le112 biomass_male_le112,cv_bio_male_le112 cv_biomass_male_le112,ci_bio_male_le112 ci_biomass_male_le112,
sum_bio_male_ge113 biomass_male_ge113,cv_bio_male_ge113 cv_biomass_male_ge113,ci_bio_male_ge113 ci_biomass_male_ge113,
sum_bio_male_ge120 biomass_male_ge120,cv_bio_male_ge120 cv_biomass_male_ge120,ci_bio_male_ge120 ci_biomass_male_ge120,
sum_bio_male_ge125 biomass_male_ge125,cv_bio_male_ge125 cv_biomass_male_ge125,ci_bio_male_ge125 ci_biomass_male_ge125,
sum_bio_male_113to124 biomass_male_113to124,cv_bio_male_113to124 cv_biomass_male_113to124,ci_bio_male_113to124 ci_biomass_male_113to124,
sum_bio_male_total biomass_male_total,cv_bio_male_total cv_biomass_male_total,ci_bio_male_total ci_biomass_male_total,
sum_bio_female_immature biomass_female_immature,cv_bio_female_immature cv_biomass_female_immature,ci_bio_female_immature ci_biomass_female_immature,
sum_bio_female_barren biomass_female_barren,cv_bio_female_barren cv_biomass_female_barren,ci_bio_female_barren ci_biomass_female_barren,
sum_bio_female_oldnoegg biomass_female_oldnoegg,cv_bio_female_oldnoegg cv_biomass_female_oldnoegg,ci_bio_female_oldnoegg ci_biomass_female_oldnoegg,
sum_bio_female_noegg biomass_female_noegg,cv_bio_female_noegg cv_biomass_female_noegg,ci_bio_female_noegg ci_biomass_female_noegg,
sum_bio_female_hatched biomass_female_hatched,cv_bio_female_hatched cv_biomass_female_hatched,ci_bio_female_hatched ci_biomass_female_hatched,
sum_bio_female_ne_tot biomass_female_ne_tot,cv_bio_female_ne_tot cv_biomass_female_ne_tot,ci_bio_female_ne_tot ci_biomass_female_ne_tot,
sum_bio_female_trace biomass_female_trace,cv_bio_female_trace cv_biomass_female_trace,ci_bio_female_trace ci_biomass_female_trace,
sum_bio_female_quarter biomass_female_quarter,cv_bio_female_quarter cv_biomass_female_quarter,ci_bio_female_quarter ci_biomass_female_quarter,
sum_bio_female_half biomass_female_half,cv_bio_female_half cv_biomass_female_half,ci_bio_female_half ci_biomass_female_half,
sum_bio_female_three_qrt biomass_female_three_qrt,cv_bio_female_three_qrt cv_biomass_female_three_qrt,ci_bio_female_three_qrt ci_biomass_female_three_qrt,
sum_bio_female_full biomass_female_full,cv_bio_female_full cv_biomass_female_full,ci_bio_female_full ci_biomass_female_full,
sum_bio_female_sc0 biomass_female_sc0,cv_bio_female_sc0 cv_biomass_female_sc0,ci_bio_female_sc0 ci_biomass_female_sc0,
sum_bio_female_sc1 biomass_female_sc1,cv_bio_female_sc1 cv_biomass_female_sc1,ci_bio_female_sc1 ci_biomass_female_sc1,
sum_bio_female_sc2 biomass_female_sc2,cv_bio_female_sc2 cv_biomass_female_sc2,ci_bio_female_sc2 ci_biomass_female_sc2,
sum_bio_female_sc3 biomass_female_sc3,cv_bio_female_sc3 cv_biomass_female_sc3,ci_bio_female_sc3 ci_biomass_female_sc3,
sum_bio_female_sc4 biomass_female_sc4,cv_bio_female_sc4 cv_biomass_female_sc4,ci_bio_female_sc4 ci_biomass_female_sc4,
sum_bio_female_sc5 biomass_female_sc5,cv_bio_female_sc5 cv_biomass_female_sc5,ci_bio_female_sc5 ci_biomass_female_sc5,
sum_bio_female_ovigerous biomass_female_ovigerous,cv_bio_female_ovigerous cv_biomass_female_ovigerous,ci_bio_female_ovigerous ci_biomass_female_ovigerous,
sum_bio_female_total biomass_female_total,cv_bio_female_total cv_biomass_female_total,ci_bio_female_total ci_biomass_female_total
from cb_bioall_sizegroup a, cb_sizegroup_confidence_bio b,cb_bio_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
order by a.survey_year;


drop table cb_ebs_pop_sc_cs_v2024;

create table cb_ebs_pop_sc_cs_v2024 as
select a.survey_year,
sum_pop_male_le112 num_male_le112,cv_pop_male_le112 cv_num_male_le112,ci_pop_male_le112 ci_num_male_le112,
sum_pop_male_ge113 num_male_ge113,cv_pop_male_ge113 cv_num_male_ge113,ci_pop_male_ge113 ci_num_male_ge113,
sum_pop_male_ge120 num_male_ge120,cv_pop_male_ge120 cv_num_male_ge120,ci_pop_male_ge120 ci_num_male_ge120,
sum_pop_male_ge125 num_male_ge125,cv_pop_male_ge125 cv_num_male_ge125,ci_pop_male_ge125 ci_num_male_ge125,
sum_pop_male_113to124 num_male_113to124,cv_pop_male_113to124 cv_num_male_113to124,ci_pop_male_113to124 ci_num_male_113to124,
sum_pop_male_total num_male_total,cv_pop_male_total cv_num_male_total,ci_pop_male_total ci_num_male_total,
sum_pop_female_immature num_female_immature,cv_pop_female_immature cv_num_female_immature,ci_pop_female_immature ci_num_female_immature,
sum_pop_female_barren num_female_barren,cv_pop_female_barren cv_num_female_barren,ci_pop_female_barren ci_num_female_barren,
sum_pop_female_oldnoegg num_female_oldnoegg,cv_pop_female_oldnoegg cv_num_female_oldnoegg,ci_pop_female_oldnoegg ci_num_female_oldnoegg,
sum_pop_female_noegg num_female_noegg,cv_pop_female_noegg cv_num_female_noegg,ci_pop_female_noegg ci_num_female_noegg,
sum_pop_female_hatched num_female_hatched,cv_pop_female_hatched cv_num_female_hatched,ci_pop_female_hatched ci_num_female_hatched,
sum_pop_female_ne_tot num_female_ne_tot,cv_pop_female_ne_tot cv_num_female_ne_tot,ci_pop_female_ne_tot ci_num_female_ne_tot,
sum_pop_female_trace num_female_trace,cv_pop_female_trace cv_num_female_trace,ci_pop_female_trace ci_num_female_trace,
sum_pop_female_quarter num_female_quarter,cv_pop_female_quarter cv_num_female_quarter,ci_pop_female_quarter ci_num_female_quarter,
sum_pop_female_half num_female_half,cv_pop_female_half cv_num_female_half,ci_pop_female_half ci_num_female_half,
sum_pop_female_three_qrt num_female_three_qrt,cv_pop_female_three_qrt cv_num_female_three_qrt,ci_pop_female_three_qrt ci_num_female_three_qrt,
sum_pop_female_full num_female_full,cv_pop_female_full cv_num_female_full,ci_pop_female_full ci_num_female_full,
sum_pop_female_sc0 num_female_sc0,cv_pop_female_sc0 cv_num_female_sc0,ci_pop_female_sc0 ci_num_female_sc0,
sum_pop_female_sc1 num_female_sc1,cv_pop_female_sc1 cv_num_female_sc1,ci_pop_female_sc1 ci_num_female_sc1,
sum_pop_female_sc2 num_female_sc2,cv_pop_female_sc2 cv_num_female_sc2,ci_pop_female_sc2 ci_num_female_sc2,
sum_pop_female_sc3 num_female_sc3,cv_pop_female_sc3 cv_num_female_sc3,ci_pop_female_sc3 ci_num_female_sc3,
sum_pop_female_sc4 num_female_sc4,cv_pop_female_sc4 cv_num_female_sc4,ci_pop_female_sc4 ci_num_female_sc4,
sum_pop_female_sc5 num_female_sc5,cv_pop_female_sc5 cv_num_female_sc5,ci_pop_female_sc5 ci_num_female_sc5,
sum_pop_female_ovigerous num_female_ovigerous,cv_pop_female_ovigerous cv_num_female_ovigerous,ci_pop_female_ovigerous ci_num_female_ovigerous,
sum_pop_female_total num_female_total,cv_pop_female_total cv_num_female_total,ci_pop_female_total ci_num_female_total
from cb_popall_sizegroup a, cb_sizegroup_confidence_pop b,cb_pop_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;


-- FOR QUADRANTS DELETE/SHUTDOWN East-166 and West-166
--  east of 166

drop table e166cb_pop_sizegroup;

create table e166cb_pop_sizegroup as
select survey_year,
sum(pop_male_le94) sum_pop_male_le94,
sum(pop_male_le102) sum_pop_male_le102,
sum(pop_male_le109) sum_pop_male_le109,
sum(pop_male_le112) sum_pop_male_le112,
sum(pop_male_le119) sum_pop_male_le119,
sum(pop_male_103to124) sum_pop_male_103to124,
sum(pop_male_113to124) sum_pop_male_113to124,
sum(pop_male_ge103) sum_pop_male_ge103,
sum(pop_male_ge113) sum_pop_male_ge113,
sum(pop_male_ge120) sum_pop_male_ge120,
sum(pop_male_ge125) sum_pop_male_ge125,
sum(pop_male_ge110) sum_pop_male_ge110,
sum(pop_male_ge138) sum_pop_male_ge138,
sum(pop_male_total) sum_pop_male_total,
sum(pop_female_immature) sum_pop_female_immature,
sum(pop_female_barren) sum_pop_female_barren,
sum(pop_female_oldnoegg) sum_pop_female_oldnoegg,
sum(pop_female_noegg) sum_pop_female_noegg,
sum(pop_female_hatched) sum_pop_female_hatched,
sum(pop_female_ne_tot) sum_pop_female_ne_tot,
sum(pop_female_trace) sum_pop_female_trace,
sum(pop_female_quarter) sum_pop_female_quarter,
sum(pop_female_half) sum_pop_female_half,
sum(pop_female_three_qrt) sum_pop_female_three_qrt,
sum(pop_female_full) sum_pop_female_full,
sum(pop_female_sc0) sum_pop_female_sc0,
sum(pop_female_sc1) sum_pop_female_sc1,
sum(pop_female_sc2) sum_pop_female_sc2,
sum(pop_female_sc3) sum_pop_female_sc3,
sum(pop_female_sc4) sum_pop_female_sc4,
sum(pop_female_sc5) sum_pop_female_sc5,
sum(pop_female_ovigerous) sum_pop_female_ovigerous,
sum(pop_female_total) sum_pop_female_total,
sum(pop_unsexed_total) sum_pop_unsexed_total,
sum(pop_gtotal) sum_pop_gtotal
from cb_popbystratum_sizegroup
where district = 'East 166'
group by survey_year
order by survey_year;

drop table e166cb_bio_sizegroup;

create table e166cb_bio_sizegroup as
select survey_year,
sum(bio_male_le94) sum_bio_male_le94,
sum(bio_male_le102) sum_bio_male_le102,
sum(bio_male_le109) sum_bio_male_le109,
sum(bio_male_le112) sum_bio_male_le112,
sum(bio_male_le119) sum_bio_male_le119,
sum(bio_male_103to124) sum_bio_male_103to124,
sum(bio_male_113to124) sum_bio_male_113to124,
sum(bio_male_ge103) sum_bio_male_ge103,
sum(bio_male_ge113) sum_bio_male_ge113,
sum(bio_male_ge120) sum_bio_male_ge120,
sum(bio_male_ge125) sum_bio_male_ge125,
sum(bio_male_ge110) sum_bio_male_ge110,
sum(bio_male_ge138) sum_bio_male_ge138,
sum(bio_male_total) sum_bio_male_total,
sum(bio_female_immature) sum_bio_female_immature,
sum(bio_female_barren) sum_bio_female_barren,
sum(bio_female_oldnoegg) sum_bio_female_oldnoegg,
sum(bio_female_noegg) sum_bio_female_noegg,
sum(bio_female_hatched) sum_bio_female_hatched,
sum(bio_female_ne_tot) sum_bio_female_ne_tot,
sum(bio_female_trace) sum_bio_female_trace,
sum(bio_female_quarter) sum_bio_female_quarter,
sum(bio_female_half) sum_bio_female_half,
sum(bio_female_three_qrt) sum_bio_female_three_qrt,
sum(bio_female_full) sum_bio_female_full,
sum(bio_female_sc0) sum_bio_female_sc0,
sum(bio_female_sc1) sum_bio_female_sc1,
sum(bio_female_sc2) sum_bio_female_sc2,
sum(bio_female_sc3) sum_bio_female_sc3,
sum(bio_female_sc4) sum_bio_female_sc4,
sum(bio_female_sc5) sum_bio_female_sc5,
sum(bio_female_ovigerous) sum_bio_female_ovigerous,
sum(bio_female_total) sum_bio_female_total,
sum(bio_gtotal) sum_bio_gtotal
from cb_biobystratum_sizegroup
where district = 'East 166'
group by survey_year
order by survey_year;


drop table e166cb_varpop_sizegroup_sum;

create table e166cb_varpop_sizegroup_sum as
select distinct survey_year,
sum(varpop_male_le94) sum_varpop_male_le94,
sum(varpop_male_le102) sum_varpop_male_le102,
sum(varpop_male_le109) sum_varpop_male_le109,
sum(varpop_male_le112) sum_varpop_male_le112,
sum(varpop_male_le119) sum_varpop_male_le119,
sum(varpop_male_103to124) sum_varpop_male_103to124,
sum(varpop_male_113to124) sum_varpop_male_113to124,
sum(varpop_male_ge103) sum_varpop_male_ge103,
sum(varpop_male_ge113) sum_varpop_male_ge113,
sum(varpop_male_ge120) sum_varpop_male_ge120,
sum(varpop_male_ge125) sum_varpop_male_ge125,
sum(varpop_male_ge110) sum_varpop_male_ge110,
sum(varpop_male_ge138) sum_varpop_male_ge138,
sum(varpop_male_total) sum_varpop_male_total,
sum(varpop_female_immature) sum_varpop_female_immature,
sum(varpop_female_barren) sum_varpop_female_barren,
sum(varpop_female_oldnoegg) sum_varpop_female_oldnoegg,
sum(varpop_female_noegg) sum_varpop_female_noegg,
sum(varpop_female_hatched) sum_varpop_female_hatched,
sum(varpop_female_ne_tot) sum_varpop_female_ne_tot,
sum(varpop_female_trace) sum_varpop_female_trace,
sum(varpop_female_quarter) sum_varpop_female_quarter,
sum(varpop_female_half) sum_varpop_female_half,
sum(varpop_female_three_qrt) sum_varpop_female_three_qrt,
sum(varpop_female_full) sum_varpop_female_full,
sum(varpop_female_sc0) sum_varpop_female_sc0,
sum(varpop_female_sc1) sum_varpop_female_sc1,
sum(varpop_female_sc2) sum_varpop_female_sc2,
sum(varpop_female_sc3) sum_varpop_female_sc3,
sum(varpop_female_sc4) sum_varpop_female_sc4,
sum(varpop_female_sc5) sum_varpop_female_sc5,
sum(varpop_female_ovigerous) sum_varpop_female_ovigerous,
sum(varpop_female_total) sum_varpop_female_total,
sum(varpop_unsexed_total) sum_varpop_unsexed_total,
sum(varpop_gtotal) sum_varpop_gtotal
from cb_variancepop_sizegroup
where district = 'East 166'
group by survey_year
order by survey_year;

drop table e166cb_varbio_sizegroup_sum;

create table e166cb_varbio_sizegroup_sum as
select distinct survey_year,
sum(varbio_male_le94) sum_varbio_male_le94,
sum(varbio_male_le102) sum_varbio_male_le102,
sum(varbio_male_le109) sum_varbio_male_le109,
sum(varbio_male_le112) sum_varbio_male_le112,
sum(varbio_male_le119) sum_varbio_male_le119,
sum(varbio_male_103to124) sum_varbio_male_103to124,
sum(varbio_male_113to124) sum_varbio_male_113to124,
sum(varbio_male_ge103) sum_varbio_male_ge103,
sum(varbio_male_ge113) sum_varbio_male_ge113,
sum(varbio_male_ge120) sum_varbio_male_ge120,
sum(varbio_male_ge125) sum_varbio_male_ge125,
sum(varbio_male_ge110) sum_varbio_male_ge110,
sum(varbio_male_ge138) sum_varbio_male_ge138,
sum(varbio_male_total) sum_varbio_male_total,
sum(varbio_female_immature) sum_varbio_female_immature,
sum(varbio_female_barren) sum_varbio_female_barren,
sum(varbio_female_oldnoegg) sum_varbio_female_oldnoegg,
sum(varbio_female_noegg) sum_varbio_female_noegg,
sum(varbio_female_hatched) sum_varbio_female_hatched,
sum(varbio_female_ne_tot) sum_varbio_female_ne_tot,
sum(varbio_female_trace) sum_varbio_female_trace,
sum(varbio_female_quarter) sum_varbio_female_quarter,
sum(varbio_female_half) sum_varbio_female_half,
sum(varbio_female_three_qrt) sum_varbio_female_three_qrt,
sum(varbio_female_full) sum_varbio_female_full,
sum(varbio_female_sc0) sum_varbio_female_sc0,
sum(varbio_female_sc1) sum_varbio_female_sc1,
sum(varbio_female_sc2) sum_varbio_female_sc2,
sum(varbio_female_sc3) sum_varbio_female_sc3,
sum(varbio_female_sc4) sum_varbio_female_sc4,
sum(varbio_female_sc5) sum_varbio_female_sc5,
sum(varbio_female_ovigerous) sum_varbio_female_ovigerous,
sum(varbio_female_total) sum_varbio_female_total,
sum(varbio_gtotal) sum_varbio_gtotal
from cb_variancebio_sizegroup
where district = 'East 166'
group by survey_year
order by survey_year;

drop table e166cb_pop_sizegroup_cv;

create table e166cb_pop_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_pop_male_le94 <> 0
	 then ((sqrt(sum_varpop_male_le94))/sum_pop_male_le94)
	 else 0
	 end) cv_pop_male_le94,
(CASE
	 when sum_pop_male_le102 <> 0
	 then ((sqrt(sum_varpop_male_le102))/sum_pop_male_le102)
	 else 0
	 end) cv_pop_male_le102,   
(CASE
	 when sum_pop_male_le109 <> 0
	 then ((sqrt(sum_varpop_male_le109))/sum_pop_male_le109)
	 else 0
	 end) cv_pop_male_le109,
(CASE
	 when sum_pop_male_le112 <> 0
	 then ((sqrt(sum_varpop_male_le112))/sum_pop_male_le112)
	 else 0
	 end) cv_pop_male_le112,
(CASE
	 when sum_pop_male_le119 <> 0
	 then ((sqrt(sum_varpop_male_le119))/sum_pop_male_le119)
	 else 0
	 end) cv_pop_male_le119,
(CASE
	 when sum_pop_male_103to124 <> 0
	 then ((sqrt(sum_varpop_male_103to124))/sum_pop_male_103to124)
	 else 0
	 end) cv_pop_male_103to124,   
(CASE
	 when sum_pop_male_113to124 <> 0
	 then ((sqrt(sum_varpop_male_113to124))/sum_pop_male_113to124)
	 else 0
	 end) cv_pop_male_113to124,
(CASE
	 when sum_pop_male_ge103 <> 0
	 then ((sqrt(sum_varpop_male_ge103))/sum_pop_male_ge103)
	 else 0
	 end) cv_pop_male_ge103,	       
(CASE
	 when sum_pop_male_ge113 <> 0
	 then ((sqrt(sum_varpop_male_ge113))/sum_pop_male_ge113)
	 else 0
	 end) cv_pop_male_ge113,	    
(CASE
	 when sum_pop_male_ge120 <> 0
	 then ((sqrt(sum_varpop_male_ge120))/sum_pop_male_ge120)
	 else 0
	 end) cv_pop_male_ge120,	 
(CASE
	 when sum_pop_male_ge125 <> 0
	 then ((sqrt(sum_varpop_male_ge125))/sum_pop_male_ge125)
	 else 0
	 end) cv_pop_male_ge125,	 	    
(CASE
	 when sum_pop_male_ge110 <> 0
	 then ((sqrt(sum_varpop_male_ge110))/sum_pop_male_ge110)
	 else 0
	 end) cv_pop_male_ge110,	 	 
(CASE
	 when sum_pop_male_ge138 <> 0
	 then ((sqrt(sum_varpop_male_ge138))/sum_pop_male_ge138)
	 else 0
	 end) cv_pop_male_ge138,
(CASE
	 when sum_pop_male_total <> 0
	 then ((sqrt(sum_varpop_male_total))/sum_pop_male_total)
	 else 0
	 end) cv_pop_male_total,
(CASE
	 when sum_pop_female_immature <> 0
	 then ((sqrt(sum_varpop_female_immature))/sum_pop_female_immature)
	 else 0
	 end) cv_pop_female_immature,	     
(CASE
	 when sum_pop_female_barren <> 0
	 then ((sqrt(sum_varpop_female_barren))/sum_pop_female_barren)
	 else 0
	 end) cv_pop_female_barren,
(CASE
	 when sum_pop_female_oldnoegg <> 0
	 then ((sqrt(sum_varpop_female_oldnoegg))/sum_pop_female_oldnoegg)
	 else 0
	 end) cv_pop_female_oldnoegg,
(CASE
	 when sum_pop_female_noegg <> 0
	 then ((sqrt(sum_varpop_female_noegg))/sum_pop_female_noegg)
	 else 0
	 end) cv_pop_female_noegg,
(CASE
	 when sum_pop_female_hatched <> 0
	 then ((sqrt(sum_varpop_female_hatched))/sum_pop_female_hatched)
	 else 0
	 end) cv_pop_female_hatched,	  
(CASE
	 when sum_pop_female_ne_tot <> 0
	 then ((sqrt(sum_varpop_female_ne_tot))/sum_pop_female_ne_tot)
	 else 0
	 end) cv_pop_female_ne_tot,	  
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
	 when sum_pop_female_three_qrt <> 0
	 then ((sqrt(sum_varpop_female_three_qrt))/sum_pop_female_three_qrt)
	 else 0
	 end) cv_pop_female_three_qrt,
(CASE
	 when sum_pop_female_full <> 0
	 then ((sqrt(sum_varpop_female_full))/sum_pop_female_full)
	 else 0
	 end) cv_pop_female_full,
(CASE
	 when sum_pop_female_sc0 <> 0
	 then ((sqrt(sum_varpop_female_sc0))/sum_pop_female_sc0)
	 else 0
	 end) cv_pop_female_sc0,
(CASE
	 when sum_pop_female_sc1 <> 0
	 then ((sqrt(sum_varpop_female_sc1))/sum_pop_female_sc1)
	 else 0
	 end) cv_pop_female_sc1,
(CASE
	 when sum_pop_female_sc2 <> 0
	 then ((sqrt(sum_varpop_female_sc2))/sum_pop_female_sc2)
	 else 0
	 end) cv_pop_female_sc2,
(CASE
	 when sum_pop_female_sc3 <> 0
	 then ((sqrt(sum_varpop_female_sc3))/sum_pop_female_sc3)
	 else 0
	 end) cv_pop_female_sc3,
(CASE
	 when sum_pop_female_sc4 <> 0
	 then ((sqrt(sum_varpop_female_sc4))/sum_pop_female_sc4)
	 else 0
	 end) cv_pop_female_sc4,
(CASE
	 when sum_pop_female_sc5 <> 0
	 then ((sqrt(sum_varpop_female_sc5))/sum_pop_female_sc5)
	 else 0
	 end) cv_pop_female_sc5,
(CASE
	 when sum_pop_female_ovigerous <> 0
	 then ((sqrt(sum_varpop_female_ovigerous))/sum_pop_female_ovigerous)
	 else 0
	 end) cv_pop_female_ovigerous,	  
(CASE
	 when sum_pop_female_total <> 0
	 then ((sqrt(sum_varpop_female_total))/sum_pop_female_total)
	 else 0
	 end) cv_pop_female_total,
(CASE
	 when sum_pop_unsexed_total <> 0
	 then ((sqrt(sum_varpop_unsexed_total))/sum_pop_unsexed_total)
	 else 0
	 end) cv_pop_unsexed_total,
(CASE
	 when sum_pop_gtotal <> 0
	 then ((sqrt(sum_varpop_gtotal))/sum_pop_gtotal)
	 else 0
	 end) cv_pop_gtotal	 	  	 	 	 	  
from e166cb_varpop_sizegroup_sum a, e166cb_pop_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;	 	 	 

drop table e166cb_bio_sizegroup_cv;

create table e166cb_bio_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_bio_male_le94 <> 0
	 then ((sqrt(sum_varbio_male_le94))/sum_bio_male_le94)
	 else 0
	 end) cv_bio_male_le94,
(CASE
	 when sum_bio_male_le102 <> 0
	 then ((sqrt(sum_varbio_male_le102))/sum_bio_male_le102)
	 else 0
	 end) cv_bio_male_le102,   
(CASE
	 when sum_bio_male_le109 <> 0
	 then ((sqrt(sum_varbio_male_le109))/sum_bio_male_le109)
	 else 0
	 end) cv_bio_male_le109,
(CASE
	 when sum_bio_male_le112 <> 0
	 then ((sqrt(sum_varbio_male_le112))/sum_bio_male_le112)
	 else 0
	 end) cv_bio_male_le112,
(CASE
	 when sum_bio_male_le119 <> 0
	 then ((sqrt(sum_varbio_male_le119))/sum_bio_male_le119)
	 else 0
	 end) cv_bio_male_le119,
(CASE
	 when sum_bio_male_103to124 <> 0
	 then ((sqrt(sum_varbio_male_103to124))/sum_bio_male_103to124)
	 else 0
	 end) cv_bio_male_103to124,
(CASE
	 when sum_bio_male_113to124 <> 0
	 then ((sqrt(sum_varbio_male_113to124))/sum_bio_male_113to124)
	 else 0
	 end) cv_bio_male_113to124,   
(CASE
	 when sum_bio_male_ge103 <> 0
	 then ((sqrt(sum_varbio_male_ge103))/sum_bio_male_ge103)
	 else 0
	 end) cv_bio_male_ge103,	       
(CASE
	 when sum_bio_male_ge113 <> 0
	 then ((sqrt(sum_varbio_male_ge113))/sum_bio_male_ge113)
	 else 0
	 end) cv_bio_male_ge113,	    
(CASE
	 when sum_bio_male_ge120 <> 0
	 then ((sqrt(sum_varbio_male_ge120))/sum_bio_male_ge120)
	 else 0
	 end) cv_bio_male_ge120,	 
(CASE
	 when sum_bio_male_ge125 <> 0
	 then ((sqrt(sum_varbio_male_ge125))/sum_bio_male_ge125)
	 else 0
	 end) cv_bio_male_ge125,	 	    
(CASE
	 when sum_bio_male_ge110 <> 0
	 then ((sqrt(sum_varbio_male_ge110))/sum_bio_male_ge110)
	 else 0
	 end) cv_bio_male_ge110,	 	 
(CASE
	 when sum_bio_male_ge138 <> 0
	 then ((sqrt(sum_varbio_male_ge138))/sum_bio_male_ge138)
	 else 0
	 end) cv_bio_male_ge138,
(CASE
	 when sum_bio_male_total <> 0
	 then ((sqrt(sum_varbio_male_total))/sum_bio_male_total)
	 else 0
	 end) cv_bio_male_total,
(CASE
	 when sum_bio_female_immature <> 0
	 then ((sqrt(sum_varbio_female_immature))/sum_bio_female_immature)
	 else 0
	 end) cv_bio_female_immature,	     
(CASE
	 when sum_bio_female_barren <> 0
	 then ((sqrt(sum_varbio_female_barren))/sum_bio_female_barren)
	 else 0
	 end) cv_bio_female_barren,
(CASE
	 when sum_bio_female_oldnoegg <> 0
	 then ((sqrt(sum_varbio_female_oldnoegg))/sum_bio_female_oldnoegg)
	 else 0
	 end) cv_bio_female_oldnoegg,
(CASE
	 when sum_bio_female_noegg <> 0
	 then ((sqrt(sum_varbio_female_noegg))/sum_bio_female_noegg)
	 else 0
	 end) cv_bio_female_noegg,
(CASE
	 when sum_bio_female_hatched <> 0
	 then ((sqrt(sum_varbio_female_hatched))/sum_bio_female_hatched)
	 else 0
	 end) cv_bio_female_hatched,
(CASE
	 when sum_bio_female_ne_tot <> 0
	 then ((sqrt(sum_varbio_female_ne_tot))/sum_bio_female_ne_tot)
	 else 0
	 end) cv_bio_female_ne_tot,
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
	 when sum_bio_female_three_qrt <> 0
	 then ((sqrt(sum_varbio_female_three_qrt))/sum_bio_female_three_qrt)
	 else 0
	 end) cv_bio_female_three_qrt,
(CASE
	 when sum_bio_female_full <> 0
	 then ((sqrt(sum_varbio_female_full))/sum_bio_female_full)
	 else 0
	 end) cv_bio_female_full,
(CASE
	 when sum_bio_female_sc0 <> 0
	 then ((sqrt(sum_varbio_female_sc0))/sum_bio_female_sc0)
	 else 0
	 end) cv_bio_female_sc0,
(CASE
	 when sum_bio_female_sc1 <> 0
	 then ((sqrt(sum_varbio_female_sc1))/sum_bio_female_sc1)
	 else 0
	 end) cv_bio_female_sc1,
(CASE
	 when sum_bio_female_sc2 <> 0
	 then ((sqrt(sum_varbio_female_sc2))/sum_bio_female_sc2)
	 else 0
	 end) cv_bio_female_sc2,
(CASE
	 when sum_bio_female_sc3 <> 0
	 then ((sqrt(sum_varbio_female_sc3))/sum_bio_female_sc3)
	 else 0
	 end) cv_bio_female_sc3,
(CASE
	 when sum_bio_female_sc4 <> 0
	 then ((sqrt(sum_varbio_female_sc4))/sum_bio_female_sc4)
	 else 0
	 end) cv_bio_female_sc4,
(CASE
	 when sum_bio_female_sc5 <> 0
	 then ((sqrt(sum_varbio_female_sc5))/sum_bio_female_sc5)
	 else 0
	 end) cv_bio_female_sc5,
(CASE
	 when sum_bio_female_ovigerous <> 0
	 then ((sqrt(sum_varbio_female_ovigerous))/sum_bio_female_ovigerous)
	 else 0
	 end) cv_bio_female_ovigerous,	
(CASE
	 when sum_bio_female_total <> 0
	 then ((sqrt(sum_varbio_female_total))/sum_bio_female_total)
	 else 0
	 end) cv_bio_female_total,
(CASE
	 when sum_bio_gtotal <> 0
	 then ((sqrt(sum_varbio_gtotal))/sum_bio_gtotal)
	 else 0
	 end) cv_bio_gtotal 
from e166cb_varbio_sizegroup_sum a, e166cb_bio_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;


-- CI calcs

create or replace view e166cb_sizegroup_stderr_pop as
select distinct survey_year,
(sqrt(sum_varpop_male_le94)) stderr_pop_male_le94,
(sqrt(sum_varpop_male_le102)) stderr_pop_male_le102,
(sqrt(sum_varpop_male_le109)) stderr_pop_male_le109,
(sqrt(sum_varpop_male_le112)) stderr_pop_male_le112,
(sqrt(sum_varpop_male_le119)) stderr_pop_male_le119,
(sqrt(sum_varpop_male_103to124)) stderr_pop_male_103to124,
(sqrt(sum_varpop_male_113to124)) stderr_pop_male_113to124,
(sqrt(sum_varpop_male_ge103)) stderr_pop_male_ge103,
(sqrt(sum_varpop_male_ge113)) stderr_pop_male_ge113,
(sqrt(sum_varpop_male_ge120)) stderr_pop_male_ge120,
(sqrt(sum_varpop_male_ge125)) stderr_pop_male_ge125,
(sqrt(sum_varpop_male_ge110)) stderr_pop_male_ge110,
(sqrt(sum_varpop_male_ge138)) stderr_pop_male_ge138,
(sqrt(sum_varpop_male_total)) stderr_pop_male_total,
(sqrt(sum_varpop_female_immature)) stderr_pop_female_immature,
(sqrt(sum_varpop_female_barren)) stderr_pop_female_barren,
(sqrt(sum_varpop_female_oldnoegg)) stderr_pop_female_oldnoegg,
(sqrt(sum_varpop_female_noegg)) stderr_pop_female_noegg,
(sqrt(sum_varpop_female_hatched)) stderr_pop_female_hatched,
(sqrt(sum_varpop_female_ne_tot)) stderr_pop_female_ne_tot,
(sqrt(sum_varpop_female_trace)) stderr_pop_female_trace,
(sqrt(sum_varpop_female_quarter)) stderr_pop_female_quarter,
(sqrt(sum_varpop_female_half)) stderr_pop_female_half,
(sqrt(sum_varpop_female_three_qrt)) stderr_pop_female_three_qrt,
(sqrt(sum_varpop_female_full)) stderr_pop_female_full,
(sqrt(sum_varpop_female_sc0)) stderr_pop_female_sc0,
(sqrt(sum_varpop_female_sc1)) stderr_pop_female_sc1,
(sqrt(sum_varpop_female_sc2)) stderr_pop_female_sc2,
(sqrt(sum_varpop_female_sc3)) stderr_pop_female_sc3,
(sqrt(sum_varpop_female_sc4)) stderr_pop_female_sc4,
(sqrt(sum_varpop_female_sc5)) stderr_pop_female_sc5,
(sqrt(sum_varpop_female_ovigerous)) stderr_pop_female_ovigerous,
(sqrt(sum_varpop_female_total)) stderr_pop_female_total,
(sqrt(sum_varpop_unsexed_total)) stderr_pop_unsexed_total,
(sqrt(sum_varpop_gtotal)) stderr_pop_gtotal
from e166cb_varpop_sizegroup_sum;

create or replace view e166cb_sizegroup_stderr_bio as
select distinct survey_year,
(sqrt(sum_varbio_male_le94)) stderr_bio_male_le94,
(sqrt(sum_varbio_male_le102)) stderr_bio_male_le102,
(sqrt(sum_varbio_male_le109)) stderr_bio_male_le109,
(sqrt(sum_varbio_male_le112)) stderr_bio_male_le112,
(sqrt(sum_varbio_male_le119)) stderr_bio_male_le119,
(sqrt(sum_varbio_male_103to124)) stderr_bio_male_103to124,
(sqrt(sum_varbio_male_113to124)) stderr_bio_male_113to124,
(sqrt(sum_varbio_male_ge103)) stderr_bio_male_ge103,
(sqrt(sum_varbio_male_ge113)) stderr_bio_male_ge113,
(sqrt(sum_varbio_male_ge120)) stderr_bio_male_ge120,
(sqrt(sum_varbio_male_ge125)) stderr_bio_male_ge125,
(sqrt(sum_varbio_male_ge110)) stderr_bio_male_ge110,
(sqrt(sum_varbio_male_ge138)) stderr_bio_male_ge138,
(sqrt(sum_varbio_male_total)) stderr_bio_male_total,
(sqrt(sum_varbio_female_immature)) stderr_bio_female_immature,
(sqrt(sum_varbio_female_barren)) stderr_bio_female_barren,
(sqrt(sum_varbio_female_oldnoegg)) stderr_bio_female_oldnoegg,
(sqrt(sum_varbio_female_noegg)) stderr_bio_female_noegg,
(sqrt(sum_varbio_female_hatched)) stderr_bio_female_hatched,
(sqrt(sum_varbio_female_ne_tot)) stderr_bio_female_ne_tot,
(sqrt(sum_varbio_female_trace)) stderr_bio_female_trace,
(sqrt(sum_varbio_female_quarter)) stderr_bio_female_quarter,
(sqrt(sum_varbio_female_half)) stderr_bio_female_half,
(sqrt(sum_varbio_female_three_qrt)) stderr_bio_female_three_qrt,
(sqrt(sum_varbio_female_full)) stderr_bio_female_full,
(sqrt(sum_varbio_female_sc0)) stderr_bio_female_sc0,
(sqrt(sum_varbio_female_sc1)) stderr_bio_female_sc1,
(sqrt(sum_varbio_female_sc2)) stderr_bio_female_sc2,
(sqrt(sum_varbio_female_sc3)) stderr_bio_female_sc3,
(sqrt(sum_varbio_female_sc4)) stderr_bio_female_sc4,
(sqrt(sum_varbio_female_sc5)) stderr_bio_female_sc5,
(sqrt(sum_varbio_female_ovigerous)) stderr_bio_female_ovigerous,
(sqrt(sum_varbio_female_total)) stderr_bio_female_total,
(sqrt(sum_varbio_gtotal)) stderr_bio_gtotal
from e166cb_varbio_sizegroup_sum;

drop table e166cb_sizegroup_ci_pop;

create table e166cb_sizegroup_ci_pop as
select distinct survey_year,
(1.96 * stderr_pop_male_le94) ci_pop_male_le94,
(1.96 * stderr_pop_male_le102) ci_pop_male_le102,
(1.96 * stderr_pop_male_le109) ci_pop_male_le109,
(1.96 * stderr_pop_male_le112) ci_pop_male_le112,
(1.96 * stderr_pop_male_le119) ci_pop_male_le119,
(1.96 * stderr_pop_male_103to124) ci_pop_male_103to124,
(1.96 * stderr_pop_male_113to124) ci_pop_male_113to124,
(1.96 * stderr_pop_male_ge103) ci_pop_male_ge103,
(1.96 * stderr_pop_male_ge113) ci_pop_male_ge113,
(1.96 * stderr_pop_male_ge120) ci_pop_male_ge120,
(1.96 * stderr_pop_male_ge125) ci_pop_male_ge125,
(1.96 * stderr_pop_male_ge110) ci_pop_male_ge110,
(1.96 * stderr_pop_male_ge138) ci_pop_male_ge138,
(1.96 * stderr_pop_male_total) ci_pop_male_total,
(1.96 * stderr_pop_female_immature) ci_pop_female_immature,
(1.96 * stderr_pop_female_barren) ci_pop_female_barren,
(1.96 * stderr_pop_female_oldnoegg) ci_pop_female_oldnoegg,
(1.96 * stderr_pop_female_noegg) ci_pop_female_noegg,
(1.96 * stderr_pop_female_hatched) ci_pop_female_hatched,
(1.96 * stderr_pop_female_ne_tot) ci_pop_female_ne_tot,
(1.96 * stderr_pop_female_trace) ci_pop_female_trace,
(1.96 * stderr_pop_female_quarter) ci_pop_female_quarter,
(1.96 * stderr_pop_female_half) ci_pop_female_half,
(1.96 * stderr_pop_female_three_qrt) ci_pop_female_three_qrt,
(1.96 * stderr_pop_female_full) ci_pop_female_full,
(1.96 * stderr_pop_female_sc0) ci_pop_female_sc0,
(1.96 * stderr_pop_female_sc1) ci_pop_female_sc1,
(1.96 * stderr_pop_female_sc2) ci_pop_female_sc2,
(1.96 * stderr_pop_female_sc3) ci_pop_female_sc3,
(1.96 * stderr_pop_female_sc4) ci_pop_female_sc4,
(1.96 * stderr_pop_female_sc5) ci_pop_female_sc5,
(1.96 * stderr_pop_female_ovigerous) ci_pop_female_ovigerous,
(1.96 * stderr_pop_female_total) ci_pop_female_total,
(1.96 * stderr_pop_unsexed_total) ci_pop_unsexed_total,
(1.96 * stderr_pop_gtotal) ci_pop_gtotal
from e166cb_sizegroup_stderr_pop;

drop table e166cb_sizegroup_ci_bio;

create table e166cb_sizegroup_ci_bio as
select distinct survey_year,
((1.96 * stderr_bio_male_le94)) ci_bio_male_le94,
((1.96 * stderr_bio_male_le102)) ci_bio_male_le102,
((1.96 * stderr_bio_male_le109)) ci_bio_male_le109,
((1.96 * stderr_bio_male_le112)) ci_bio_male_le112,
((1.96 * stderr_bio_male_le119)) ci_bio_male_le119,
((1.96 * stderr_bio_male_103to124)) ci_bio_male_103to124,
((1.96 * stderr_bio_male_113to124)) ci_bio_male_113to124,
((1.96 * stderr_bio_male_ge103)) ci_bio_male_ge103,
((1.96 * stderr_bio_male_ge113)) ci_bio_male_ge113,
((1.96 * stderr_bio_male_ge120)) ci_bio_male_ge120,
((1.96 * stderr_bio_male_ge125)) ci_bio_male_ge125,
((1.96 * stderr_bio_male_ge110)) ci_bio_male_ge110,
((1.96 * stderr_bio_male_ge138)) ci_bio_male_ge138,
((1.96 * stderr_bio_male_total)) ci_bio_male_total,
((1.96 * stderr_bio_female_immature)) ci_bio_female_immature,
((1.96 * stderr_bio_female_barren)) ci_bio_female_barren,
((1.96 * stderr_bio_female_oldnoegg)) ci_bio_female_oldnoegg,
((1.96 * stderr_bio_female_noegg)) ci_bio_female_noegg,
((1.96 * stderr_bio_female_hatched)) ci_bio_female_hatched,
((1.96 * stderr_bio_female_ne_tot)) ci_bio_female_ne_tot,
((1.96 * stderr_bio_female_trace)) ci_bio_female_trace,
((1.96 * stderr_bio_female_quarter)) ci_bio_female_quarter,
((1.96 * stderr_bio_female_half)) ci_bio_female_half,
((1.96 * stderr_bio_female_three_qrt)) ci_bio_female_three_qrt,
((1.96 * stderr_bio_female_full)) ci_bio_female_full,
((1.96 * stderr_bio_female_sc0)) ci_bio_female_sc0,
((1.96 * stderr_bio_female_sc1)) ci_bio_female_sc1,
((1.96 * stderr_bio_female_sc2)) ci_bio_female_sc2,
((1.96 * stderr_bio_female_sc3)) ci_bio_female_sc3,
((1.96 * stderr_bio_female_sc4)) ci_bio_female_sc4,
((1.96 * stderr_bio_female_sc5)) ci_bio_female_sc5,
((1.96 * stderr_bio_female_ovigerous)) ci_bio_female_ovigerous,
((1.96 * stderr_bio_female_total)) ci_bio_female_total,
((1.96 * stderr_bio_gtotal)) ci_bio_gtotal
from e166cb_sizegroup_stderr_bio;


--  west of 166W

drop table w166cb_pop_sizegroup;

create table w166cb_pop_sizegroup as
select survey_year,
sum(pop_male_le94) sum_pop_male_le94,
sum(pop_male_le102) sum_pop_male_le102,
sum(pop_male_le109) sum_pop_male_le109,
sum(pop_male_le112) sum_pop_male_le112,
sum(pop_male_le119) sum_pop_male_le119,
sum(pop_male_103to124) sum_pop_male_103to124,
sum(pop_male_113to124) sum_pop_male_113to124,
sum(pop_male_ge103) sum_pop_male_ge103,
sum(pop_male_ge113) sum_pop_male_ge113,
sum(pop_male_ge120) sum_pop_male_ge120,
sum(pop_male_ge125) sum_pop_male_ge125,
sum(pop_male_ge110) sum_pop_male_ge110,
sum(pop_male_ge138) sum_pop_male_ge138,
sum(pop_male_total) sum_pop_male_total,
sum(pop_female_immature) sum_pop_female_immature,
sum(pop_female_barren) sum_pop_female_barren,
sum(pop_female_oldnoegg) sum_pop_female_oldnoegg,
sum(pop_female_noegg) sum_pop_female_noegg,
sum(pop_female_hatched) sum_pop_female_hatched,
sum(pop_female_ne_tot) sum_pop_female_ne_tot,
sum(pop_female_trace) sum_pop_female_trace,
sum(pop_female_quarter) sum_pop_female_quarter,
sum(pop_female_half) sum_pop_female_half,
sum(pop_female_three_qrt) sum_pop_female_three_qrt,
sum(pop_female_full) sum_pop_female_full,
sum(pop_female_sc0) sum_pop_female_sc0,
sum(pop_female_sc1) sum_pop_female_sc1,
sum(pop_female_sc2) sum_pop_female_sc2,
sum(pop_female_sc3) sum_pop_female_sc3,
sum(pop_female_sc4) sum_pop_female_sc4,
sum(pop_female_sc5) sum_pop_female_sc5,
sum(pop_female_ovigerous) sum_pop_female_ovigerous,
sum(pop_female_total) sum_pop_female_total,
sum(pop_unsexed_total) sum_pop_unsexed_total,
sum(pop_gtotal) sum_pop_gtotal
from cb_popbystratum_sizegroup
where district in ('Pribilof MTCA','St. Matthew MTCA','West 166')
group by survey_year
order by survey_year;

drop table w166cb_bio_sizegroup;

create table w166cb_bio_sizegroup as
select survey_year,
sum(bio_male_le94) sum_bio_male_le94,
sum(bio_male_le102) sum_bio_male_le102,
sum(bio_male_le109) sum_bio_male_le109,
sum(bio_male_le112) sum_bio_male_le112,
sum(bio_male_le119) sum_bio_male_le119,
sum(bio_male_103to124) sum_bio_male_103to124,
sum(bio_male_113to124) sum_bio_male_113to124,
sum(bio_male_ge103) sum_bio_male_ge103,
sum(bio_male_ge113) sum_bio_male_ge113,
sum(bio_male_ge120) sum_bio_male_ge120,
sum(bio_male_ge125) sum_bio_male_ge125,
sum(bio_male_ge110) sum_bio_male_ge110,
sum(bio_male_ge138) sum_bio_male_ge138,
sum(bio_male_total) sum_bio_male_total,
sum(bio_female_immature) sum_bio_female_immature,
sum(bio_female_barren) sum_bio_female_barren,
sum(bio_female_oldnoegg) sum_bio_female_oldnoegg,
sum(bio_female_noegg) sum_bio_female_noegg,
sum(bio_female_hatched) sum_bio_female_hatched,
sum(bio_female_ne_tot) sum_bio_female_ne_tot,
sum(bio_female_trace) sum_bio_female_trace,
sum(bio_female_quarter) sum_bio_female_quarter,
sum(bio_female_half) sum_bio_female_half,
sum(bio_female_three_qrt) sum_bio_female_three_qrt,
sum(bio_female_full) sum_bio_female_full,
sum(bio_female_sc0) sum_bio_female_sc0,
sum(bio_female_sc1) sum_bio_female_sc1,
sum(bio_female_sc2) sum_bio_female_sc2,
sum(bio_female_sc3) sum_bio_female_sc3,
sum(bio_female_sc4) sum_bio_female_sc4,
sum(bio_female_sc5) sum_bio_female_sc5,
sum(bio_female_ovigerous) sum_bio_female_ovigerous,
sum(bio_female_total) sum_bio_female_total,
sum(bio_gtotal) sum_bio_gtotal
from cb_biobystratum_sizegroup
where district in ('Pribilof MTCA','St. Matthew MTCA','West 166')
group by survey_year
order by survey_year;


drop table w166cb_varpop_sizegroup_sum;

create table w166cb_varpop_sizegroup_sum as
select distinct survey_year,
sum(varpop_male_le94) sum_varpop_male_le94,
sum(varpop_male_le102) sum_varpop_male_le102,
sum(varpop_male_le109) sum_varpop_male_le109,
sum(varpop_male_le112) sum_varpop_male_le112,
sum(varpop_male_le119) sum_varpop_male_le119,
sum(varpop_male_103to124) sum_varpop_male_103to124,
sum(varpop_male_113to124) sum_varpop_male_113to124,
sum(varpop_male_ge103) sum_varpop_male_ge103,
sum(varpop_male_ge113) sum_varpop_male_ge113,
sum(varpop_male_ge120) sum_varpop_male_ge120,
sum(varpop_male_ge125) sum_varpop_male_ge125,
sum(varpop_male_ge110) sum_varpop_male_ge110,
sum(varpop_male_ge138) sum_varpop_male_ge138,
sum(varpop_male_total) sum_varpop_male_total,
sum(varpop_female_immature) sum_varpop_female_immature,
sum(varpop_female_barren) sum_varpop_female_barren,
sum(varpop_female_oldnoegg) sum_varpop_female_oldnoegg,
sum(varpop_female_noegg) sum_varpop_female_noegg,
sum(varpop_female_hatched) sum_varpop_female_hatched,
sum(varpop_female_ne_tot) sum_varpop_female_ne_tot,
sum(varpop_female_trace) sum_varpop_female_trace,
sum(varpop_female_quarter) sum_varpop_female_quarter,
sum(varpop_female_half) sum_varpop_female_half,
sum(varpop_female_three_qrt) sum_varpop_female_three_qrt,
sum(varpop_female_full) sum_varpop_female_full,
sum(varpop_female_sc0) sum_varpop_female_sc0,
sum(varpop_female_sc1) sum_varpop_female_sc1,
sum(varpop_female_sc2) sum_varpop_female_sc2,
sum(varpop_female_sc3) sum_varpop_female_sc3,
sum(varpop_female_sc4) sum_varpop_female_sc4,
sum(varpop_female_sc5) sum_varpop_female_sc5,
sum(varpop_female_ovigerous) sum_varpop_female_ovigerous,
sum(varpop_female_total) sum_varpop_female_total,
sum(varpop_unsexed_total) sum_varpop_unsexed_total,
sum(varpop_gtotal) sum_varpop_gtotal
from cb_variancepop_sizegroup
where district in ('Pribilof MTCA','St. Matthew MTCA','West 166')
group by survey_year
order by survey_year;

drop table w166cb_varbio_sizegroup_sum;

create table w166cb_varbio_sizegroup_sum as
select distinct survey_year,
sum(varbio_male_le94) sum_varbio_male_le94,
sum(varbio_male_le102) sum_varbio_male_le102,
sum(varbio_male_le109) sum_varbio_male_le109,
sum(varbio_male_le112) sum_varbio_male_le112,
sum(varbio_male_le119) sum_varbio_male_le119,
sum(varbio_male_103to124) sum_varbio_male_103to124,
sum(varbio_male_113to124) sum_varbio_male_113to124,
sum(varbio_male_ge103) sum_varbio_male_ge103,
sum(varbio_male_ge113) sum_varbio_male_ge113,
sum(varbio_male_ge120) sum_varbio_male_ge120,
sum(varbio_male_ge125) sum_varbio_male_ge125,
sum(varbio_male_ge110) sum_varbio_male_ge110,
sum(varbio_male_ge138) sum_varbio_male_ge138,
sum(varbio_male_total) sum_varbio_male_total,
sum(varbio_female_immature) sum_varbio_female_immature,
sum(varbio_female_barren) sum_varbio_female_barren,
sum(varbio_female_oldnoegg) sum_varbio_female_oldnoegg,
sum(varbio_female_noegg) sum_varbio_female_noegg,
sum(varbio_female_hatched) sum_varbio_female_hatched,
sum(varbio_female_ne_tot) sum_varbio_female_ne_tot,
sum(varbio_female_trace) sum_varbio_female_trace,
sum(varbio_female_quarter) sum_varbio_female_quarter,
sum(varbio_female_half) sum_varbio_female_half,
sum(varbio_female_three_qrt) sum_varbio_female_three_qrt,
sum(varbio_female_full) sum_varbio_female_full,
sum(varbio_female_sc0) sum_varbio_female_sc0,
sum(varbio_female_sc1) sum_varbio_female_sc1,
sum(varbio_female_sc2) sum_varbio_female_sc2,
sum(varbio_female_sc3) sum_varbio_female_sc3,
sum(varbio_female_sc4) sum_varbio_female_sc4,
sum(varbio_female_sc5) sum_varbio_female_sc5,
sum(varbio_female_ovigerous) sum_varbio_female_ovigerous,
sum(varbio_female_total) sum_varbio_female_total,
sum(varbio_gtotal) sum_varbio_gtotal
from cb_variancebio_sizegroup
where district in ('Pribilof MTCA','St. Matthew MTCA','West 166')
group by survey_year
order by survey_year;

drop table w166cb_pop_sizegroup_cv;

create table w166cb_pop_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_pop_male_le94 <> 0
	 then ((sqrt(sum_varpop_male_le94))/sum_pop_male_le94)
	 else 0
	 end) cv_pop_male_le94,
(CASE
	 when sum_pop_male_le102 <> 0
	 then ((sqrt(sum_varpop_male_le102))/sum_pop_male_le102)
	 else 0
	 end) cv_pop_male_le102,   
(CASE
	 when sum_pop_male_le109 <> 0
	 then ((sqrt(sum_varpop_male_le109))/sum_pop_male_le109)
	 else 0
	 end) cv_pop_male_le109,
(CASE
	 when sum_pop_male_le112 <> 0
	 then ((sqrt(sum_varpop_male_le112))/sum_pop_male_le112)
	 else 0
	 end) cv_pop_male_le112,
(CASE
	 when sum_pop_male_le119 <> 0
	 then ((sqrt(sum_varpop_male_le119))/sum_pop_male_le119)
	 else 0
	 end) cv_pop_male_le119,
(CASE
	 when sum_pop_male_103to124 <> 0
	 then ((sqrt(sum_varpop_male_103to124))/sum_pop_male_103to124)
	 else 0
	 end) cv_pop_male_103to124,   
(CASE
	 when sum_pop_male_113to124 <> 0
	 then ((sqrt(sum_varpop_male_113to124))/sum_pop_male_113to124)
	 else 0
	 end) cv_pop_male_113to124,
(CASE
	 when sum_pop_male_ge103 <> 0
	 then ((sqrt(sum_varpop_male_ge103))/sum_pop_male_ge103)
	 else 0
	 end) cv_pop_male_ge103,	       
(CASE
	 when sum_pop_male_ge113 <> 0
	 then ((sqrt(sum_varpop_male_ge113))/sum_pop_male_ge113)
	 else 0
	 end) cv_pop_male_ge113,	    
(CASE
	 when sum_pop_male_ge120 <> 0
	 then ((sqrt(sum_varpop_male_ge120))/sum_pop_male_ge120)
	 else 0
	 end) cv_pop_male_ge120,	 
(CASE
	 when sum_pop_male_ge125 <> 0
	 then ((sqrt(sum_varpop_male_ge125))/sum_pop_male_ge125)
	 else 0
	 end) cv_pop_male_ge125,	 	    
(CASE
	 when sum_pop_male_ge110 <> 0
	 then ((sqrt(sum_varpop_male_ge110))/sum_pop_male_ge110)
	 else 0
	 end) cv_pop_male_ge110,	 	 
(CASE
	 when sum_pop_male_ge138 <> 0
	 then ((sqrt(sum_varpop_male_ge138))/sum_pop_male_ge138)
	 else 0
	 end) cv_pop_male_ge138,
(CASE
	 when sum_pop_male_total <> 0
	 then ((sqrt(sum_varpop_male_total))/sum_pop_male_total)
	 else 0
	 end) cv_pop_male_total,
(CASE
	 when sum_pop_female_immature <> 0
	 then ((sqrt(sum_varpop_female_immature))/sum_pop_female_immature)
	 else 0
	 end) cv_pop_female_immature,	     
(CASE
	 when sum_pop_female_barren <> 0
	 then ((sqrt(sum_varpop_female_barren))/sum_pop_female_barren)
	 else 0
	 end) cv_pop_female_barren,
(CASE
	 when sum_pop_female_oldnoegg <> 0
	 then ((sqrt(sum_varpop_female_oldnoegg))/sum_pop_female_oldnoegg)
	 else 0
	 end) cv_pop_female_oldnoegg,
(CASE
	 when sum_pop_female_noegg <> 0
	 then ((sqrt(sum_varpop_female_noegg))/sum_pop_female_noegg)
	 else 0
	 end) cv_pop_female_noegg,
(CASE
	 when sum_pop_female_hatched <> 0
	 then ((sqrt(sum_varpop_female_hatched))/sum_pop_female_hatched)
	 else 0
	 end) cv_pop_female_hatched,	  
(CASE
	 when sum_pop_female_ne_tot <> 0
	 then ((sqrt(sum_varpop_female_ne_tot))/sum_pop_female_ne_tot)
	 else 0
	 end) cv_pop_female_ne_tot,	  
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
	 when sum_pop_female_three_qrt <> 0
	 then ((sqrt(sum_varpop_female_three_qrt))/sum_pop_female_three_qrt)
	 else 0
	 end) cv_pop_female_three_qrt,
(CASE
	 when sum_pop_female_full <> 0
	 then ((sqrt(sum_varpop_female_full))/sum_pop_female_full)
	 else 0
	 end) cv_pop_female_full,
(CASE
	 when sum_pop_female_sc0 <> 0
	 then ((sqrt(sum_varpop_female_sc0))/sum_pop_female_sc0)
	 else 0
	 end) cv_pop_female_sc0,
(CASE
	 when sum_pop_female_sc1 <> 0
	 then ((sqrt(sum_varpop_female_sc1))/sum_pop_female_sc1)
	 else 0
	 end) cv_pop_female_sc1,
(CASE
	 when sum_pop_female_sc2 <> 0
	 then ((sqrt(sum_varpop_female_sc2))/sum_pop_female_sc2)
	 else 0
	 end) cv_pop_female_sc2,
(CASE
	 when sum_pop_female_sc3 <> 0
	 then ((sqrt(sum_varpop_female_sc3))/sum_pop_female_sc3)
	 else 0
	 end) cv_pop_female_sc3,
(CASE
	 when sum_pop_female_sc4 <> 0
	 then ((sqrt(sum_varpop_female_sc4))/sum_pop_female_sc4)
	 else 0
	 end) cv_pop_female_sc4,
(CASE
	 when sum_pop_female_sc5 <> 0
	 then ((sqrt(sum_varpop_female_sc5))/sum_pop_female_sc5)
	 else 0
	 end) cv_pop_female_sc5,
(CASE
	 when sum_pop_female_ovigerous <> 0
	 then ((sqrt(sum_varpop_female_ovigerous))/sum_pop_female_ovigerous)
	 else 0
	 end) cv_pop_female_ovigerous,	  
(CASE
	 when sum_pop_female_total <> 0
	 then ((sqrt(sum_varpop_female_total))/sum_pop_female_total)
	 else 0
	 end) cv_pop_female_total,
(CASE
	 when sum_pop_unsexed_total <> 0
	 then ((sqrt(sum_varpop_unsexed_total))/sum_pop_unsexed_total)
	 else 0
	 end) cv_pop_unsexed_total,
(CASE
	 when sum_pop_gtotal <> 0
	 then ((sqrt(sum_varpop_gtotal))/sum_pop_gtotal)
	 else 0
	 end) cv_pop_gtotal	  	 	 	 	  	 	 	 	  
from w166cb_varpop_sizegroup_sum a, w166cb_pop_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;	 	 	 

drop table w166cb_bio_sizegroup_cv;

create table w166cb_bio_sizegroup_cv as
select a.survey_year,
(CASE
	 when sum_bio_male_le94 <> 0
	 then ((sqrt(sum_varbio_male_le94))/sum_bio_male_le94)
	 else 0
	 end) cv_bio_male_le94,
(CASE
	 when sum_bio_male_le102 <> 0
	 then ((sqrt(sum_varbio_male_le102))/sum_bio_male_le102)
	 else 0
	 end) cv_bio_male_le102,   
(CASE
	 when sum_bio_male_le109 <> 0
	 then ((sqrt(sum_varbio_male_le109))/sum_bio_male_le109)
	 else 0
	 end) cv_bio_male_le109,
(CASE
	 when sum_bio_male_le112 <> 0
	 then ((sqrt(sum_varbio_male_le112))/sum_bio_male_le112)
	 else 0
	 end) cv_bio_male_le112,
(CASE
	 when sum_bio_male_le119 <> 0
	 then ((sqrt(sum_varbio_male_le119))/sum_bio_male_le119)
	 else 0
	 end) cv_bio_male_le119,
(CASE
	 when sum_bio_male_103to124 <> 0
	 then ((sqrt(sum_varbio_male_103to124))/sum_bio_male_103to124)
	 else 0
	 end) cv_bio_male_103to124,
(CASE
	 when sum_bio_male_113to124 <> 0
	 then ((sqrt(sum_varbio_male_113to124))/sum_bio_male_113to124)
	 else 0
	 end) cv_bio_male_113to124,   
(CASE
	 when sum_bio_male_ge103 <> 0
	 then ((sqrt(sum_varbio_male_ge103))/sum_bio_male_ge103)
	 else 0
	 end) cv_bio_male_ge103,	       
(CASE
	 when sum_bio_male_ge113 <> 0
	 then ((sqrt(sum_varbio_male_ge113))/sum_bio_male_ge113)
	 else 0
	 end) cv_bio_male_ge113,	    
(CASE
	 when sum_bio_male_ge120 <> 0
	 then ((sqrt(sum_varbio_male_ge120))/sum_bio_male_ge120)
	 else 0
	 end) cv_bio_male_ge120,	 
(CASE
	 when sum_bio_male_ge125 <> 0
	 then ((sqrt(sum_varbio_male_ge125))/sum_bio_male_ge125)
	 else 0
	 end) cv_bio_male_ge125,	 	    
(CASE
	 when sum_bio_male_ge110 <> 0
	 then ((sqrt(sum_varbio_male_ge110))/sum_bio_male_ge110)
	 else 0
	 end) cv_bio_male_ge110,	 	 
(CASE
	 when sum_bio_male_ge138 <> 0
	 then ((sqrt(sum_varbio_male_ge138))/sum_bio_male_ge138)
	 else 0
	 end) cv_bio_male_ge138,
(CASE
	 when sum_bio_male_total <> 0
	 then ((sqrt(sum_varbio_male_total))/sum_bio_male_total)
	 else 0
	 end) cv_bio_male_total,
(CASE
	 when sum_bio_female_immature <> 0
	 then ((sqrt(sum_varbio_female_immature))/sum_bio_female_immature)
	 else 0
	 end) cv_bio_female_immature,	     
(CASE
	 when sum_bio_female_barren <> 0
	 then ((sqrt(sum_varbio_female_barren))/sum_bio_female_barren)
	 else 0
	 end) cv_bio_female_barren,
(CASE
	 when sum_bio_female_oldnoegg <> 0
	 then ((sqrt(sum_varbio_female_oldnoegg))/sum_bio_female_oldnoegg)
	 else 0
	 end) cv_bio_female_oldnoegg,
(CASE
	 when sum_bio_female_noegg <> 0
	 then ((sqrt(sum_varbio_female_noegg))/sum_bio_female_noegg)
	 else 0
	 end) cv_bio_female_noegg,
(CASE
	 when sum_bio_female_hatched <> 0
	 then ((sqrt(sum_varbio_female_hatched))/sum_bio_female_hatched)
	 else 0
	 end) cv_bio_female_hatched,
(CASE
	 when sum_bio_female_ne_tot <> 0
	 then ((sqrt(sum_varbio_female_ne_tot))/sum_bio_female_ne_tot)
	 else 0
	 end) cv_bio_female_ne_tot,
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
	 when sum_bio_female_three_qrt <> 0
	 then ((sqrt(sum_varbio_female_three_qrt))/sum_bio_female_three_qrt)
	 else 0
	 end) cv_bio_female_three_qrt,
(CASE
	 when sum_bio_female_full <> 0
	 then ((sqrt(sum_varbio_female_full))/sum_bio_female_full)
	 else 0
	 end) cv_bio_female_full,
(CASE
	 when sum_bio_female_sc0 <> 0
	 then ((sqrt(sum_varbio_female_sc0))/sum_bio_female_sc0)
	 else 0
	 end) cv_bio_female_sc0,
(CASE
	 when sum_bio_female_sc1 <> 0
	 then ((sqrt(sum_varbio_female_sc1))/sum_bio_female_sc1)
	 else 0
	 end) cv_bio_female_sc1,
(CASE
	 when sum_bio_female_sc2 <> 0
	 then ((sqrt(sum_varbio_female_sc2))/sum_bio_female_sc2)
	 else 0
	 end) cv_bio_female_sc2,
(CASE
	 when sum_bio_female_sc3 <> 0
	 then ((sqrt(sum_varbio_female_sc3))/sum_bio_female_sc3)
	 else 0
	 end) cv_bio_female_sc3,
(CASE
	 when sum_bio_female_sc4 <> 0
	 then ((sqrt(sum_varbio_female_sc4))/sum_bio_female_sc4)
	 else 0
	 end) cv_bio_female_sc4,
(CASE
	 when sum_bio_female_sc5 <> 0
	 then ((sqrt(sum_varbio_female_sc5))/sum_bio_female_sc5)
	 else 0
	 end) cv_bio_female_sc5,
(CASE
	 when sum_bio_female_ovigerous <> 0
	 then ((sqrt(sum_varbio_female_ovigerous))/sum_bio_female_ovigerous)
	 else 0
	 end) cv_bio_female_ovigerous,	
(CASE
	 when sum_bio_female_total <> 0
	 then ((sqrt(sum_varbio_female_total))/sum_bio_female_total)
	 else 0
	 end) cv_bio_female_total,
(CASE
	 when sum_bio_gtotal <> 0
	 then ((sqrt(sum_varbio_gtotal))/sum_bio_gtotal)
	 else 0
	 end) cv_bio_gtotal	  	
from w166cb_varbio_sizegroup_sum a,w166cb_bio_sizegroup b
where a.survey_year = b.survey_year
order by a.survey_year;


-- CI calcs

create or replace view w166cb_sizegroup_stderr_pop as
select distinct survey_year,
(sqrt(sum_varpop_male_le94)) stderr_pop_male_le94,
(sqrt(sum_varpop_male_le102)) stderr_pop_male_le102,
(sqrt(sum_varpop_male_le109)) stderr_pop_male_le109,
(sqrt(sum_varpop_male_le112)) stderr_pop_male_le112,
(sqrt(sum_varpop_male_le119)) stderr_pop_male_le119,
(sqrt(sum_varpop_male_103to124)) stderr_pop_male_103to124,
(sqrt(sum_varpop_male_113to124)) stderr_pop_male_113to124,
(sqrt(sum_varpop_male_ge103)) stderr_pop_male_ge103,
(sqrt(sum_varpop_male_ge113)) stderr_pop_male_ge113,
(sqrt(sum_varpop_male_ge120)) stderr_pop_male_ge120,
(sqrt(sum_varpop_male_ge125)) stderr_pop_male_ge125,
(sqrt(sum_varpop_male_ge110)) stderr_pop_male_ge110,
(sqrt(sum_varpop_male_ge138)) stderr_pop_male_ge138,
(sqrt(sum_varpop_male_total)) stderr_pop_male_total,
(sqrt(sum_varpop_female_immature)) stderr_pop_female_immature,
(sqrt(sum_varpop_female_barren)) stderr_pop_female_barren,
(sqrt(sum_varpop_female_oldnoegg)) stderr_pop_female_oldnoegg,
(sqrt(sum_varpop_female_noegg)) stderr_pop_female_noegg,
(sqrt(sum_varpop_female_hatched)) stderr_pop_female_hatched,
(sqrt(sum_varpop_female_ne_tot)) stderr_pop_female_ne_tot,
(sqrt(sum_varpop_female_trace)) stderr_pop_female_trace,
(sqrt(sum_varpop_female_quarter)) stderr_pop_female_quarter,
(sqrt(sum_varpop_female_half)) stderr_pop_female_half,
(sqrt(sum_varpop_female_three_qrt)) stderr_pop_female_three_qrt,
(sqrt(sum_varpop_female_full)) stderr_pop_female_full,
(sqrt(sum_varpop_female_sc0)) stderr_pop_female_sc0,
(sqrt(sum_varpop_female_sc1)) stderr_pop_female_sc1,
(sqrt(sum_varpop_female_sc2)) stderr_pop_female_sc2,
(sqrt(sum_varpop_female_sc3)) stderr_pop_female_sc3,
(sqrt(sum_varpop_female_sc4)) stderr_pop_female_sc4,
(sqrt(sum_varpop_female_sc5)) stderr_pop_female_sc5,
(sqrt(sum_varpop_female_ovigerous)) stderr_pop_female_ovigerous,
(sqrt(sum_varpop_female_total)) stderr_pop_female_total,
(sqrt(sum_varpop_unsexed_total)) stderr_pop_unsexed_total,
(sqrt(sum_varpop_gtotal)) stderr_pop_gtotal
from w166cb_varpop_sizegroup_sum;

create or replace view w166cb_sizegroup_stderr_bio as
select distinct survey_year,
(sqrt(sum_varbio_male_le94)) stderr_bio_male_le94,
(sqrt(sum_varbio_male_le102)) stderr_bio_male_le102,
(sqrt(sum_varbio_male_le109)) stderr_bio_male_le109,
(sqrt(sum_varbio_male_le112)) stderr_bio_male_le112,
(sqrt(sum_varbio_male_le119)) stderr_bio_male_le119,
(sqrt(sum_varbio_male_103to124)) stderr_bio_male_103to124,
(sqrt(sum_varbio_male_113to124)) stderr_bio_male_113to124,
(sqrt(sum_varbio_male_ge103)) stderr_bio_male_ge103,
(sqrt(sum_varbio_male_ge113)) stderr_bio_male_ge113,
(sqrt(sum_varbio_male_ge120)) stderr_bio_male_ge120,
(sqrt(sum_varbio_male_ge125)) stderr_bio_male_ge125,
(sqrt(sum_varbio_male_ge110)) stderr_bio_male_ge110,
(sqrt(sum_varbio_male_ge138)) stderr_bio_male_ge138,
(sqrt(sum_varbio_male_total)) stderr_bio_male_total,
(sqrt(sum_varbio_female_immature)) stderr_bio_female_immature,
(sqrt(sum_varbio_female_barren)) stderr_bio_female_barren,
(sqrt(sum_varbio_female_oldnoegg)) stderr_bio_female_oldnoegg,
(sqrt(sum_varbio_female_noegg)) stderr_bio_female_noegg,
(sqrt(sum_varbio_female_hatched)) stderr_bio_female_hatched,
(sqrt(sum_varbio_female_ne_tot)) stderr_bio_female_ne_tot,
(sqrt(sum_varbio_female_trace)) stderr_bio_female_trace,
(sqrt(sum_varbio_female_quarter)) stderr_bio_female_quarter,
(sqrt(sum_varbio_female_half)) stderr_bio_female_half,
(sqrt(sum_varbio_female_three_qrt)) stderr_bio_female_three_qrt,
(sqrt(sum_varbio_female_full)) stderr_bio_female_full,
(sqrt(sum_varbio_female_sc0)) stderr_bio_female_sc0,
(sqrt(sum_varbio_female_sc1)) stderr_bio_female_sc1,
(sqrt(sum_varbio_female_sc2)) stderr_bio_female_sc2,
(sqrt(sum_varbio_female_sc3)) stderr_bio_female_sc3,
(sqrt(sum_varbio_female_sc4)) stderr_bio_female_sc4,
(sqrt(sum_varbio_female_sc5)) stderr_bio_female_sc5,
(sqrt(sum_varbio_female_ovigerous)) stderr_bio_female_ovigerous,
(sqrt(sum_varbio_female_total)) stderr_bio_female_total,
(sqrt(sum_varbio_gtotal)) stderr_bio_gtotal
from w166cb_varbio_sizegroup_sum;

drop table w166cb_sizegroup_ci_pop;

create table w166cb_sizegroup_ci_pop as
select distinct survey_year,
(1.96 * stderr_pop_male_le94) ci_pop_male_le94,
(1.96 * stderr_pop_male_le102) ci_pop_male_le102,
(1.96 * stderr_pop_male_le109) ci_pop_male_le109,
(1.96 * stderr_pop_male_le112) ci_pop_male_le112,
(1.96 * stderr_pop_male_le119) ci_pop_male_le119,
(1.96 * stderr_pop_male_103to124) ci_pop_male_103to124,
(1.96 * stderr_pop_male_113to124) ci_pop_male_113to124,
(1.96 * stderr_pop_male_ge103) ci_pop_male_ge103,
(1.96 * stderr_pop_male_ge113) ci_pop_male_ge113,
(1.96 * stderr_pop_male_ge120) ci_pop_male_ge120,
(1.96 * stderr_pop_male_ge125) ci_pop_male_ge125,
(1.96 * stderr_pop_male_ge110) ci_pop_male_ge110,
(1.96 * stderr_pop_male_ge138) ci_pop_male_ge138,
(1.96 * stderr_pop_male_total) ci_pop_male_total,
(1.96 * stderr_pop_female_immature) ci_pop_female_immature,
(1.96 * stderr_pop_female_barren) ci_pop_female_barren,
(1.96 * stderr_pop_female_oldnoegg) ci_pop_female_oldnoegg,
(1.96 * stderr_pop_female_noegg) ci_pop_female_noegg,
(1.96 * stderr_pop_female_hatched) ci_pop_female_hatched,
(1.96 * stderr_pop_female_ne_tot) ci_pop_female_ne_tot,
(1.96 * stderr_pop_female_trace) ci_pop_female_trace,
(1.96 * stderr_pop_female_quarter) ci_pop_female_quarter,
(1.96 * stderr_pop_female_half) ci_pop_female_half,
(1.96 * stderr_pop_female_three_qrt) ci_pop_female_three_qrt,
(1.96 * stderr_pop_female_full) ci_pop_female_full,
(1.96 * stderr_pop_female_sc0) ci_pop_female_sc0,
(1.96 * stderr_pop_female_sc1) ci_pop_female_sc1,
(1.96 * stderr_pop_female_sc2) ci_pop_female_sc2,
(1.96 * stderr_pop_female_sc3) ci_pop_female_sc3,
(1.96 * stderr_pop_female_sc4) ci_pop_female_sc4,
(1.96 * stderr_pop_female_sc5) ci_pop_female_sc5,
(1.96 * stderr_pop_female_ovigerous) ci_pop_female_ovigerous,
(1.96 * stderr_pop_female_total) ci_pop_female_total,
(1.96 * stderr_pop_unsexed_total) ci_pop_unsexed_total,
(1.96 * stderr_pop_gtotal) ci_pop_gtotal
from w166cb_sizegroup_stderr_pop;

drop table w166cb_sizegroup_ci_bio;

create table w166cb_sizegroup_ci_bio as
select distinct survey_year,
((1.96 * stderr_bio_male_le94)) ci_bio_male_le94,
((1.96 * stderr_bio_male_le102)) ci_bio_male_le102,
((1.96 * stderr_bio_male_le109)) ci_bio_male_le109,
((1.96 * stderr_bio_male_le112)) ci_bio_male_le112,
((1.96 * stderr_bio_male_le119)) ci_bio_male_le119,
((1.96 * stderr_bio_male_103to124)) ci_bio_male_103to124,
((1.96 * stderr_bio_male_113to124)) ci_bio_male_113to124,
((1.96 * stderr_bio_male_ge103)) ci_bio_male_ge103,
((1.96 * stderr_bio_male_ge113)) ci_bio_male_ge113,
((1.96 * stderr_bio_male_ge120)) ci_bio_male_ge120,
((1.96 * stderr_bio_male_ge125)) ci_bio_male_ge125,
((1.96 * stderr_bio_male_ge110)) ci_bio_male_ge110,
((1.96 * stderr_bio_male_ge138)) ci_bio_male_ge138,
((1.96 * stderr_bio_male_total)) ci_bio_male_total,
((1.96 * stderr_bio_female_immature)) ci_bio_female_immature,
((1.96 * stderr_bio_female_barren)) ci_bio_female_barren,
((1.96 * stderr_bio_female_oldnoegg)) ci_bio_female_oldnoegg,
((1.96 * stderr_bio_female_noegg)) ci_bio_female_noegg,
((1.96 * stderr_bio_female_hatched)) ci_bio_female_hatched,
((1.96 * stderr_bio_female_ne_tot)) ci_bio_female_ne_tot,
((1.96 * stderr_bio_female_trace)) ci_bio_female_trace,
((1.96 * stderr_bio_female_quarter)) ci_bio_female_quarter,
((1.96 * stderr_bio_female_half)) ci_bio_female_half,
((1.96 * stderr_bio_female_three_qrt)) ci_bio_female_three_qrt,
((1.96 * stderr_bio_female_full)) ci_bio_female_full,
((1.96 * stderr_bio_female_sc0)) ci_bio_female_sc0,
((1.96 * stderr_bio_female_sc1)) ci_bio_female_sc1,
((1.96 * stderr_bio_female_sc2)) ci_bio_female_sc2,
((1.96 * stderr_bio_female_sc3)) ci_bio_female_sc3,
((1.96 * stderr_bio_female_sc4)) ci_bio_female_sc4,
((1.96 * stderr_bio_female_sc5)) ci_bio_female_sc5,
((1.96 * stderr_bio_female_ovigerous)) ci_bio_female_ovigerous,
((1.96 * stderr_bio_female_total)) ci_bio_female_total,
((1.96 * stderr_bio_gtotal)) ci_bio_gtotal
from w166cb_sizegroup_stderr_bio;


-- East of 166W

drop table cb_e166_bio_matfem_sc_cs_v2024;

create table cb_e166_bio_matfem_sc_cs_v2024 as
select a.survey_year,
sum_bio_male_le112 biomass_male_le112,cv_bio_male_le112 cv_biomass_male_le112,ci_bio_male_le112 ci_biomass_male_le112,
sum_bio_male_ge113 biomass_male_ge113,cv_bio_male_ge113 cv_biomass_male_ge113,ci_bio_male_ge113 ci_biomass_male_ge113,
sum_bio_male_ge120 biomass_male_ge120,cv_bio_male_ge120 cv_biomass_male_ge120,ci_bio_male_ge120 ci_biomass_male_ge120,
sum_bio_male_ge125 biomass_male_ge125,cv_bio_male_ge125 cv_biomass_male_ge125,ci_bio_male_ge125 ci_biomass_male_ge125,
sum_bio_male_113to124 biomass_male_113to124,cv_bio_male_113to124 cv_biomass_male_113to124,ci_bio_male_113to124 ci_biomass_male_113to124,
sum_bio_male_total biomass_male_total,cv_bio_male_total cv_biomass_male_total,ci_bio_male_total ci_biomass_male_total,
sum_bio_female_immature biomass_female_immature,cv_bio_female_immature cv_biomass_female_immature,ci_bio_female_immature ci_biomass_female_immature,
sum_bio_female_barren biomass_female_barren,cv_bio_female_barren cv_biomass_female_barren,ci_bio_female_barren ci_biomass_female_barren,
sum_bio_female_oldnoegg biomass_female_oldnoegg,cv_bio_female_oldnoegg cv_biomass_female_oldnoegg,ci_bio_female_oldnoegg ci_biomass_female_oldnoegg,
sum_bio_female_noegg biomass_female_noegg,cv_bio_female_noegg cv_biomass_female_noegg,ci_bio_female_noegg ci_biomass_female_noegg,
sum_bio_female_hatched biomass_female_hatched,cv_bio_female_hatched cv_biomass_female_hatched,ci_bio_female_hatched ci_biomass_female_hatched,
sum_bio_female_ne_tot biomass_female_ne_tot,cv_bio_female_ne_tot cv_biomass_female_ne_tot,ci_bio_female_ne_tot ci_biomass_female_ne_tot,
sum_bio_female_trace biomass_female_trace,cv_bio_female_trace cv_biomass_female_trace,ci_bio_female_trace ci_biomass_female_trace,
sum_bio_female_quarter biomass_female_quarter,cv_bio_female_quarter cv_biomass_female_quarter,ci_bio_female_quarter ci_biomass_female_quarter,
sum_bio_female_half biomass_female_half,cv_bio_female_half cv_biomass_female_half,ci_bio_female_half ci_biomass_female_half,
sum_bio_female_three_qrt biomass_female_three_qrt,cv_bio_female_three_qrt cv_biomass_female_three_qrt,ci_bio_female_three_qrt ci_biomass_female_three_qrt,
sum_bio_female_full biomass_female_full,cv_bio_female_full cv_biomass_female_full,ci_bio_female_full ci_biomass_female_full,
sum_bio_female_sc0 biomass_female_sc0,cv_bio_female_sc0 cv_biomass_female_sc0,ci_bio_female_sc0 ci_biomass_female_sc0,
sum_bio_female_sc1 biomass_female_sc1,cv_bio_female_sc1 cv_biomass_female_sc1,ci_bio_female_sc1 ci_biomass_female_sc1,
sum_bio_female_sc2 biomass_female_sc2,cv_bio_female_sc2 cv_biomass_female_sc2,ci_bio_female_sc2 ci_biomass_female_sc2,
sum_bio_female_sc3 biomass_female_sc3,cv_bio_female_sc3 cv_biomass_female_sc3,ci_bio_female_sc3 ci_biomass_female_sc3,
sum_bio_female_sc4 biomass_female_sc4,cv_bio_female_sc4 cv_biomass_female_sc4,ci_bio_female_sc4 ci_biomass_female_sc4,
sum_bio_female_sc5 biomass_female_sc5,cv_bio_female_sc5 cv_biomass_female_sc5,ci_bio_female_sc5 ci_biomass_female_sc5,
sum_bio_female_ovigerous biomass_female_ovigerous,cv_bio_female_ovigerous cv_biomass_female_ovigerous,ci_bio_female_ovigerous ci_biomass_female_ovigerous,
sum_bio_female_total biomass_female_total,cv_bio_female_total cv_biomass_female_total,ci_bio_female_total ci_biomass_female_total
from e166cb_bio_sizegroup a, e166cb_sizegroup_ci_bio b,e166cb_bio_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;


drop table cb_e166_pop_matfem_sc_cs_v2024;

create table cb_e166_pop_matfem_sc_cs_v2024 as
select a.survey_year,
sum_pop_male_le112 num_male_le112,cv_pop_male_le112 cv_num_male_le112,ci_pop_male_le112 ci_num_male_le112,
sum_pop_male_ge113 num_male_ge113,cv_pop_male_ge113 cv_num_male_ge113,ci_pop_male_ge113 ci_num_male_ge113,
sum_pop_male_ge120 num_male_ge120,cv_pop_male_ge120 cv_num_male_ge120,ci_pop_male_ge120 ci_num_male_ge120,
sum_pop_male_ge125 num_male_ge125,cv_pop_male_ge125 cv_num_male_ge125,ci_pop_male_ge125 ci_num_male_ge125,
sum_pop_male_113to124 num_male_113to124,cv_pop_male_113to124 cv_num_male_113to124,ci_pop_male_113to124 ci_num_male_113to124,
sum_pop_male_total num_male_total,cv_pop_male_total cv_num_male_total,ci_pop_male_total ci_num_male_total,
sum_pop_female_immature num_female_immature,cv_pop_female_immature cv_num_female_immature,ci_pop_female_immature ci_num_female_immature,
sum_pop_female_barren num_female_barren,cv_pop_female_barren cv_num_female_barren,ci_pop_female_barren ci_num_female_barren,
sum_pop_female_oldnoegg num_female_oldnoegg,cv_pop_female_oldnoegg cv_num_female_oldnoegg,ci_pop_female_oldnoegg ci_num_female_oldnoegg,
sum_pop_female_noegg num_female_noegg,cv_pop_female_noegg cv_num_female_noegg,ci_pop_female_noegg ci_num_female_noegg,
sum_pop_female_hatched num_female_hatched,cv_pop_female_hatched cv_num_female_hatched,ci_pop_female_hatched ci_num_female_hatched,
sum_pop_female_ne_tot num_female_ne_tot,cv_pop_female_ne_tot cv_num_female_ne_tot,ci_pop_female_ne_tot ci_num_female_ne_tot,
sum_pop_female_trace num_female_trace,cv_pop_female_trace cv_num_female_trace,ci_pop_female_trace ci_num_female_trace,
sum_pop_female_quarter num_female_quarter,cv_pop_female_quarter cv_num_female_quarter,ci_pop_female_quarter ci_num_female_quarter,
sum_pop_female_half num_female_half,cv_pop_female_half cv_num_female_half,ci_pop_female_half ci_num_female_half,
sum_pop_female_three_qrt num_female_three_qrt,cv_pop_female_three_qrt cv_num_female_three_qrt,ci_pop_female_three_qrt ci_num_female_three_qrt,
sum_pop_female_full num_female_full,cv_pop_female_full cv_num_female_full,ci_pop_female_full ci_num_female_full,
sum_pop_female_sc0 num_female_sc0,cv_pop_female_sc0 cv_num_female_sc0,ci_pop_female_sc0 ci_num_female_sc0,
sum_pop_female_sc1 num_female_sc1,cv_pop_female_sc1 cv_num_female_sc1,ci_pop_female_sc1 ci_num_female_sc1,
sum_pop_female_sc2 num_female_sc2,cv_pop_female_sc2 cv_num_female_sc2,ci_pop_female_sc2 ci_num_female_sc2,
sum_pop_female_sc3 num_female_sc3,cv_pop_female_sc3 cv_num_female_sc3,ci_pop_female_sc3 ci_num_female_sc3,
sum_pop_female_sc4 num_female_sc4,cv_pop_female_sc4 cv_num_female_sc4,ci_pop_female_sc4 ci_num_female_sc4,
sum_pop_female_sc5 num_female_sc5,cv_pop_female_sc5 cv_num_female_sc5,ci_pop_female_sc5 ci_num_female_sc5,
sum_pop_female_ovigerous num_female_ovigerous,cv_pop_female_ovigerous cv_num_female_ovigerous,ci_pop_female_ovigerous ci_num_female_ovigerous,
sum_pop_female_total num_female_total,cv_pop_female_total cv_num_female_total,ci_pop_female_total ci_num_female_total
from e166cb_pop_sizegroup a, e166cb_sizegroup_ci_pop b,e166cb_pop_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;


-- West of 166W

drop table cb_w166_bio_matfem_sc_cs_v2024;

create table cb_w166_bio_matfem_sc_cs_v2024 as
select a.survey_year,
sum_bio_male_le112 biomass_male_le112,cv_bio_male_le112 cv_biomass_male_le112,ci_bio_male_le112 ci_biomass_male_le112,
sum_bio_male_ge113 biomass_male_ge113,cv_bio_male_ge113 cv_biomass_male_ge113,ci_bio_male_ge113 ci_biomass_male_ge113,
sum_bio_male_ge120 biomass_male_ge120,cv_bio_male_ge120 cv_biomass_male_ge120,ci_bio_male_ge120 ci_biomass_male_ge120,
sum_bio_male_ge125 biomass_male_ge125,cv_bio_male_ge125 cv_biomass_male_ge125,ci_bio_male_ge125 ci_biomass_male_ge125,
sum_bio_male_113to124 biomass_male_113to124,cv_bio_male_113to124 cv_biomass_male_113to124,ci_bio_male_113to124 ci_biomass_male_113to124,
sum_bio_male_total biomass_male_total,cv_bio_male_total cv_biomass_male_total,ci_bio_male_total ci_biomass_male_total,
sum_bio_female_immature biomass_female_immature,cv_bio_female_immature cv_biomass_female_immature,ci_bio_female_immature ci_biomass_female_immature,
sum_bio_female_barren biomass_female_barren,cv_bio_female_barren cv_biomass_female_barren,ci_bio_female_barren ci_biomass_female_barren,
sum_bio_female_oldnoegg biomass_female_oldnoegg,cv_bio_female_oldnoegg cv_biomass_female_oldnoegg,ci_bio_female_oldnoegg ci_biomass_female_oldnoegg,
sum_bio_female_noegg biomass_female_noegg,cv_bio_female_noegg cv_biomass_female_noegg,ci_bio_female_noegg ci_biomass_female_noegg,
sum_bio_female_hatched biomass_female_hatched,cv_bio_female_hatched cv_biomass_female_hatched,ci_bio_female_hatched ci_biomass_female_hatched,
sum_bio_female_ne_tot biomass_female_ne_tot,cv_bio_female_ne_tot cv_biomass_female_ne_tot,ci_bio_female_ne_tot ci_biomass_female_ne_tot,
sum_bio_female_trace biomass_female_trace,cv_bio_female_trace cv_biomass_female_trace,ci_bio_female_trace ci_biomass_female_trace,
sum_bio_female_quarter biomass_female_quarter,cv_bio_female_quarter cv_biomass_female_quarter,ci_bio_female_quarter ci_biomass_female_quarter,
sum_bio_female_half biomass_female_half,cv_bio_female_half cv_biomass_female_half,ci_bio_female_half ci_biomass_female_half,
sum_bio_female_three_qrt biomass_female_three_qrt,cv_bio_female_three_qrt cv_biomass_female_three_qrt,ci_bio_female_three_qrt ci_biomass_female_three_qrt,
sum_bio_female_full biomass_female_full,cv_bio_female_full cv_biomass_female_full,ci_bio_female_full ci_biomass_female_full,
sum_bio_female_sc0 biomass_female_sc0,cv_bio_female_sc0 cv_biomass_female_sc0,ci_bio_female_sc0 ci_biomass_female_sc0,
sum_bio_female_sc1 biomass_female_sc1,cv_bio_female_sc1 cv_biomass_female_sc1,ci_bio_female_sc1 ci_biomass_female_sc1,
sum_bio_female_sc2 biomass_female_sc2,cv_bio_female_sc2 cv_biomass_female_sc2,ci_bio_female_sc2 ci_biomass_female_sc2,
sum_bio_female_sc3 biomass_female_sc3,cv_bio_female_sc3 cv_biomass_female_sc3,ci_bio_female_sc3 ci_biomass_female_sc3,
sum_bio_female_sc4 biomass_female_sc4,cv_bio_female_sc4 cv_biomass_female_sc4,ci_bio_female_sc4 ci_biomass_female_sc4,
sum_bio_female_sc5 biomass_female_sc5,cv_bio_female_sc5 cv_biomass_female_sc5,ci_bio_female_sc5 ci_biomass_female_sc5,
sum_bio_female_ovigerous biomass_female_ovigerous,cv_bio_female_ovigerous cv_biomass_female_ovigerous,ci_bio_female_ovigerous ci_biomass_female_ovigerous,
sum_bio_female_total biomass_female_total,cv_bio_female_total cv_biomass_female_total,ci_bio_female_total ci_biomass_female_total
from w166cb_bio_sizegroup a, w166cb_sizegroup_ci_bio b, w166cb_bio_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;


drop table cb_w166_pop_matfem_sc_cs_v2024;

create table cb_w166_pop_matfem_sc_cs_v2024 as
select a.survey_year,
sum_pop_male_le112 num_male_le112,cv_pop_male_le112 cv_num_male_le112,ci_pop_male_le112 ci_num_male_le112,
sum_pop_male_ge113 num_male_ge113,cv_pop_male_ge113 cv_num_male_ge113,ci_pop_male_ge113 ci_num_male_ge113,
sum_pop_male_ge120 num_male_ge120,cv_pop_male_ge120 cv_num_male_ge120,ci_pop_male_ge120 ci_num_male_ge120,
sum_pop_male_ge125 num_male_ge125,cv_pop_male_ge125 cv_num_male_ge125,ci_pop_male_ge125 ci_num_male_ge125,
sum_pop_male_113to124 num_male_113to124,cv_pop_male_113to124 cv_num_male_113to124,ci_pop_male_113to124 ci_num_male_113to124,
sum_pop_male_total num_male_total,cv_pop_male_total cv_num_male_total,ci_pop_male_total ci_num_male_total,
sum_pop_female_immature num_female_immature,cv_pop_female_immature cv_num_female_immature,ci_pop_female_immature ci_num_female_immature,
sum_pop_female_barren num_female_barren,cv_pop_female_barren cv_num_female_barren,ci_pop_female_barren ci_num_female_barren,
sum_pop_female_oldnoegg num_female_oldnoegg,cv_pop_female_oldnoegg cv_num_female_oldnoegg,ci_pop_female_oldnoegg ci_num_female_oldnoegg,
sum_pop_female_noegg num_female_noegg,cv_pop_female_noegg cv_num_female_noegg,ci_pop_female_noegg ci_num_female_noegg,
sum_pop_female_hatched num_female_hatched,cv_pop_female_hatched cv_num_female_hatched,ci_pop_female_hatched ci_num_female_hatched,
sum_pop_female_ne_tot num_female_ne_tot,cv_pop_female_ne_tot cv_num_female_ne_tot,ci_pop_female_ne_tot ci_num_female_ne_tot,
sum_pop_female_trace num_female_trace,cv_pop_female_trace cv_num_female_trace,ci_pop_female_trace ci_num_female_trace,
sum_pop_female_quarter num_female_quarter,cv_pop_female_quarter cv_num_female_quarter,ci_pop_female_quarter ci_num_female_quarter,
sum_pop_female_half num_female_half,cv_pop_female_half cv_num_female_half,ci_pop_female_half ci_num_female_half,
sum_pop_female_three_qrt num_female_three_qrt,cv_pop_female_three_qrt cv_num_female_three_qrt,ci_pop_female_three_qrt ci_num_female_three_qrt,
sum_pop_female_full num_female_full,cv_pop_female_full cv_num_female_full,ci_pop_female_full ci_num_female_full,
sum_pop_female_sc0 num_female_sc0,cv_pop_female_sc0 cv_num_female_sc0,ci_pop_female_sc0 ci_num_female_sc0,
sum_pop_female_sc1 num_female_sc1,cv_pop_female_sc1 cv_num_female_sc1,ci_pop_female_sc1 ci_num_female_sc1,
sum_pop_female_sc2 num_female_sc2,cv_pop_female_sc2 cv_num_female_sc2,ci_pop_female_sc2 ci_num_female_sc2,
sum_pop_female_sc3 num_female_sc3,cv_pop_female_sc3 cv_num_female_sc3,ci_pop_female_sc3 ci_num_female_sc3,
sum_pop_female_sc4 num_female_sc4,cv_pop_female_sc4 cv_num_female_sc4,ci_pop_female_sc4 ci_num_female_sc4,
sum_pop_female_sc5 num_female_sc5,cv_pop_female_sc5 cv_num_female_sc5,ci_pop_female_sc5 ci_num_female_sc5,
sum_pop_female_ovigerous num_female_ovigerous,cv_pop_female_ovigerous cv_num_female_ovigerous,ci_pop_female_ovigerous ci_num_female_ovigerous,
sum_pop_female_total num_female_total,cv_pop_female_total cv_num_female_total,ci_pop_female_total ci_num_female_total
from w166cb_pop_sizegroup a, w166cb_sizegroup_ci_pop b, w166cb_pop_sizegroup_cv c
where a.survey_year = b.survey_year
and a.survey_year = c.survey_year
--and a.survey_year = 2013
order by a.survey_year;



