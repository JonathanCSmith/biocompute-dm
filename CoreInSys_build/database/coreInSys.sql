


drop database coreInSys;


create database coreInSys;

use coreInSys;

create table masterProject
(
masterProjectID int unsigned auto_increment,
status varchar(20) not null,
projectName varchar(40),
projectLead varchar(40),
description text,
openDate date,
lastUpdate date,
primary key (masterProjectID)
);

create table typeLinker
(
linkID int unsigned auto_increment,
type enum('masterProject','sequencing','flow cytometry','analysis'),
parentID int unsigned,
childID int unsigned,
primary key (linkID)
);

create table state
(
stateID int unsigned auto_increment,
itemID int unsigned,
state enum('RED','GREEN','BLUE','AMBER','BROWN'),
type enum('masterProject','sequencing','flow cytometry','analysis'),
primary key (stateID)
);


create table flowCellProject
(
flowCellProjectID int unsigned auto_increment,
flowCellProjectName varchar(40),
primary key (flowCellProjectID)
);



create table seqProject
(
seqProjectID int unsigned auto_increment,
seqProjectName varchar(40),
masterProjectID int unsigned,
seqExptID varchar(10),
customerID int,
primary key (seqProjectID)
);


create table seqExperiment
(
seqExptID int unsigned auto_increment,
seqExptName varchar(40),
flowcellID varchar(20),
startDate date,
completionDate date,
genomicsLead varchar(20),
dataLocation varchar(40),
indexTagCycles tinyint,
readCycles tinyint,
pairedEnd enum('Y','N'),
primary key (seqExptID)
);

create table customer
(
customerID int,
name varchar(40),
email varchar(50),
tel varchar(30),
primary key (customerID)
);

create table laneData
(
laneID int unsigned auto_increment,
laneNumber tinyint,
seqProjectID int unsigned,
sequencingConc float,
read1ClusterDensity int,
PhiXspiked float,
spike varchar(20),
spikeRatio float,
primary key (laneID)
);

create table sampleData
(
sampleID int unsigned auto_increment,
sampleName varchar(50),
tagID int unsigned,
laneID int unsigned,
tagSequence varchar(20),
tagKit varchar(50),
analysisID int unsigned,
adaptorSequence varchar(200),
primary key (sampleID)

);

create table IDtagLibs
(
tagID int unsigned auto_increment,
libraryName varchar(40),
libraryTagID varchar(20),
tagSequence varchar(40),
primary key (tagID)

);

create table demultiplex
(
demuxID int unsigned auto_increment,
seqProjectID int unsigned,
status enum('setup','running','complete','finished'),
location varchar(500),
sourceLocation varchar(500),
JID int unsigned,
primary key (demuxID)
);






