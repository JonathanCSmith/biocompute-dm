drop database fastQCDB;


create database fastQCDB;

use fastQCDB;

create table sample
(
fastqFileID int unsigned auto_increment,
sampleID int unsigned,
primary key (fastqFileID)
);


create table summary
(
fastqFileID int unsigned not null,
Basic_Statistics enum('PASS','WARN','FAIL'),
Per_base_sequence_quality enum('PASS','WARN','FAIL'),
Per_sequence_quality_scores enum('PASS','WARN','FAIL'),
Per_base_sequence_content enum('PASS','WARN','FAIL'),
Per_base_GC_content enum('PASS','WARN','FAIL'),
Per_base_N_content enum('PASS','WARN','FAIL'),
Sequence_Length_Distribution enum('PASS','WARN','FAIL'),
Sequence_Duplication_Levels enum('PASS','WARN','FAIL'),
Overrepresented_sequences enum('PASS','WARN','FAIL'),
Kmer_Content enum('PASS','WARN','FAIL'),
primary key (fastqFileID)
);


create table Per_base_sequence_quality
(
fastqFileID int unsigned not null,
BaseStart tinyint,
BaseEnd tinyint,
Mean double,
Median float(2,1),
Lower_Quartile float(2,1),
Upper_Quartile float(2,1),
10th_Percentile float(2,1),
90th_Percentile float(2,1),
primary key (fastqFileID)
);


create table Per_sequence_quality_scores
(
fastqFileID int unsigned not null,
Quality tinyint,
Count int unsigned,
primary key (fastqFileID)
);

create table Per_base_sequence_content
(
fastqFileID int unsigned not null,
BaseStart tinyint,
BaseEnd tinyint,
G double,
A double,
T double,
C double,
primary key (fastqFileID)
);

create table Per_base_GC_content
(
fastqFileID int unsigned not null,
BaseStart tinyint,
BaseEnd tinyint,
PercGC double,
primary key (fastqFileID)
);

create table Per_sequence_GC_content
(
fastqFileID int unsigned not null,
GC_Content tinyint,
Count int unsigned,
primary key (fastqFileID)
);

create table Per_base_N_content
(
fastqFileID int unsigned not null,
BaseStart tinyint,
BaseEnd tinyint,
N_Count double,
primary key (fastqFileID)
);


create table Sequence_Length_Distribution
(
fastqFileID int unsigned not null,
Length tinyint,
Count double,
primary key (fastqFileID)
);


create table Sequence_Duplication_Levels
(
fastqFileID int unsigned not null,
Duplication_Level tinyint,
Relative_count double,
Total_Duplicate_Percentage double,
primary key (fastqFileID)
);

create table Kmer_Content
(
fastqFileID int unsigned not null,
Sequence varchar(10),
Count int unsigned,
Obs_over_Exp_Overall float,
Obs_over_Exp_Max float,
Max_Obs_over_Expected_PositionStart tinyint,
Max_Obs_over_Expected_PositionEnd tinyint,
primary key (fastqFileID)
);









