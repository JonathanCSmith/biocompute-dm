


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


create table flowCytProject
(
flowCytProjectID int unsigned auto_increment,
flowCytProjectName varchar(40),
primary key (flowCytProjectID)
);



create table seqProject
(
seqProjectID int unsigned auto_increment,
seqProjectName varchar(40),
masterProjectID int unsigned,
seqRunID varchar(10),
customerID int,
primary key (seqProjectID)
);


create table seqRun
(
seqRunID int unsigned auto_increment,
seqRunName varchar(40),
flowcellID varchar(20),
startDate date,
completionDate date,
genomicsLead varchar(20),
dataLocation varchar(40),
indexTagCycles tinyint,
readCycles tinyint,
pairedEnd enum('Y','N'),
primary key (seqRunID)
);

create table customer
(
customerID int unsigned auto_increment,
name varchar(40),
primary key (customerID)
);

create table customerContact
(
contactID int unsigned auto_increment,
customerID int unsigned,
name varchar(40),
email varchar(50),
tel varchar(30),
primary key (contactID)
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

create table fastQC
(
fastQCID int unsigned auto_increment,
status enum('setup','running','complete','finished'),
seqProjectID int unsigned,
sampleID int unsigned,
location varchar(500),
sourceLocation varchar(500),
JID int unsigned,
primary key (fastQCID)
);


create table masterProjectDocuments
(
documentID int unsigned auto_increment,
docDescription text,
docLocation varchar(500),
masterProjectID int unsigned,
primary key (documentID)
);


create table transfer
(
transferID int unsigned auto_increment,
seqProjectID int unsigned,
transLocation varchar(500),
sftpAccountID int unsigned,
dataStatus enum('origin','sftp'),
primary key (transferID)
);


create table sftpAccount
(
sftpAccountID int unsigned auto_increment,
accountLocation varchar(500),
username varchar(20),
userContact varchar(100),
creationDate date,
accountStatus enum('created','deleted','restore'),
primary key (sftpAccountID)
);





