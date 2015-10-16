-- MySQL dump 10.13  Distrib 5.1.67, for redhat-linux-gnu (x86_64)
--
-- Host: localhost    Database: coreInSys
-- ------------------------------------------------------
-- Server version	5.1.67

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `IDtagLibs`
--

DROP TABLE IF EXISTS `IDtagLibs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `IDtagLibs` (
  `tagID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `libraryName` varchar(40) DEFAULT NULL,
  `libraryTagID` varchar(20) DEFAULT NULL,
  `tagSequence` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`tagID`)
) ENGINE=MyISAM AUTO_INCREMENT=397 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `customer`
--

DROP TABLE IF EXISTS `customer`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `customer` (
  `customerID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`customerID`)
) ENGINE=MyISAM AUTO_INCREMENT=10 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `customerContact`
--

DROP TABLE IF EXISTS `customerContact`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `customerContact` (
  `contactID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `customerID` int(10) unsigned DEFAULT NULL,
  `name` varchar(40) DEFAULT NULL,
  `email` varchar(50) DEFAULT NULL,
  `tel` varchar(30) DEFAULT NULL,
  PRIMARY KEY (`contactID`)
) ENGINE=MyISAM AUTO_INCREMENT=3 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `demultiplex`
--

DROP TABLE IF EXISTS `demultiplex`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `demultiplex` (
  `demuxID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seqProjectID` int(10) unsigned DEFAULT NULL,
  `status` enum('setup','running','complete','finished') DEFAULT NULL,
  `location` varchar(500) DEFAULT NULL,
  `sourceLocation` varchar(500) DEFAULT NULL,
  `JID` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`demuxID`)
) ENGINE=MyISAM AUTO_INCREMENT=547 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `fastQC`
--

DROP TABLE IF EXISTS `fastQC`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `fastQC` (
  `fastQCID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `status` enum('setup','running','complete','finished') DEFAULT NULL,
  `seqProjectID` int(10) unsigned DEFAULT NULL,
  `sampleID` int(10) unsigned DEFAULT NULL,
  `location` varchar(500) DEFAULT NULL,
  `sourceLocation` varchar(500) DEFAULT NULL,
  `JID` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`fastQCID`)
) ENGINE=MyISAM AUTO_INCREMENT=5931 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `flowCytProject`
--

DROP TABLE IF EXISTS `flowCytProject`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `flowCytProject` (
  `flowCytProjectID` int(10) unsigned NOT NULL DEFAULT '0',
  `flowCytProjectName` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`flowCytProjectID`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `laneData`
--

DROP TABLE IF EXISTS `laneData`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `laneData` (
  `laneID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `laneNumber` tinyint(4) DEFAULT NULL,
  `seqProjectID` int(10) unsigned DEFAULT NULL,
  `sequencingConc` float DEFAULT NULL,
  `read1ClusterDensity` int(11) DEFAULT NULL,
  `PhiXspiked` float DEFAULT NULL,
  `spike` varchar(20) DEFAULT NULL,
  `spikeRatio` float DEFAULT NULL,
  PRIMARY KEY (`laneID`)
) ENGINE=MyISAM AUTO_INCREMENT=854 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `masterProject`
--

DROP TABLE IF EXISTS `masterProject`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `masterProject` (
  `masterProjectID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `status` varchar(20) NOT NULL,
  `projectName` varchar(40) DEFAULT NULL,
  `projectLead` varchar(40) DEFAULT NULL,
  `description` text,
  `openDate` date DEFAULT NULL,
  `lastUpdate` date DEFAULT NULL,
  PRIMARY KEY (`masterProjectID`)
) ENGINE=MyISAM AUTO_INCREMENT=193 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `masterProjectDocuments`
--

DROP TABLE IF EXISTS `masterProjectDocuments`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `masterProjectDocuments` (
  `documentID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `docDescription` text,
  `docLocation` varchar(500) DEFAULT NULL,
  `masterProjectID` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`documentID`)
) ENGINE=MyISAM AUTO_INCREMENT=18 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `postAlignQC`
--

DROP TABLE IF EXISTS `postAlignQC`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `postAlignQC` (
  `postAlignQCID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `status` enum('setup','running','complete','finished') DEFAULT NULL,
  `seqProjectID` int(10) unsigned DEFAULT NULL,
  `sampleID` int(10) unsigned DEFAULT NULL,
  `location` varchar(500) DEFAULT NULL,
  `sourceLocation` varchar(500) DEFAULT NULL,
  `JID` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`postAlignQCID`)
) ENGINE=MyISAM AUTO_INCREMENT=71 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sampleData`
--

DROP TABLE IF EXISTS `sampleData`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sampleData` (
  `sampleID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `sampleName` varchar(50) DEFAULT NULL,
  `tagID` int(10) unsigned DEFAULT NULL,
  `laneID` int(10) unsigned DEFAULT NULL,
  `tagSequence` varchar(20) DEFAULT NULL,
  `tagKit` varchar(50) DEFAULT NULL,
  `analysisID` int(10) unsigned DEFAULT NULL,
  `adaptorSequence` varchar(200) DEFAULT NULL,
  PRIMARY KEY (`sampleID`)
) ENGINE=MyISAM AUTO_INCREMENT=8953 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `seqProject`
--

DROP TABLE IF EXISTS `seqProject`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `seqProject` (
  `seqProjectID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seqProjectName` varchar(40) DEFAULT NULL,
  `masterProjectID` int(10) unsigned DEFAULT NULL,
  `seqRunID` int(10) unsigned DEFAULT NULL,
  `customerID` int(11) DEFAULT NULL,
  `exptType` enum('exome','RNAseq','ChIPseq','other','WGS') DEFAULT NULL,
  PRIMARY KEY (`seqProjectID`)
) ENGINE=MyISAM AUTO_INCREMENT=273 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `seqRun`
--

DROP TABLE IF EXISTS `seqRun`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `seqRun` (
  `seqRunID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seqRunName` varchar(40) DEFAULT NULL,
  `flowcellID` varchar(20) DEFAULT NULL,
  `startDate` date DEFAULT NULL,
  `completionDate` date DEFAULT NULL,
  `genomicsLead` varchar(20) DEFAULT NULL,
  `dataLocation` varchar(40) DEFAULT NULL,
  `indexTagCycles` tinyint(4) DEFAULT NULL,
  `readCycles` smallint(6) DEFAULT NULL,
  `pairedEnd` enum('Y','N') DEFAULT NULL,
  PRIMARY KEY (`seqRunID`)
) ENGINE=MyISAM AUTO_INCREMENT=137 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sftpAccount`
--

DROP TABLE IF EXISTS `sftpAccount`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sftpAccount` (
  `sftpAccountID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `accountLocation` varchar(500) DEFAULT NULL,
  `username` varchar(20) DEFAULT NULL,
  `userContact` varchar(100) DEFAULT NULL,
  `creationDate` date DEFAULT NULL,
  `accountStatus` enum('created','deleted','restore') DEFAULT NULL,
  PRIMARY KEY (`sftpAccountID`)
) ENGINE=MyISAM AUTO_INCREMENT=251 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `state`
--

DROP TABLE IF EXISTS `state`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `state` (
  `stateID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `itemID` int(10) unsigned DEFAULT NULL,
  `state` enum('RED','GREEN','BLUE','AMBER','BROWN') DEFAULT NULL,
  `type` enum('masterProject','sequencing','flow cytometry','analysis') DEFAULT NULL,
  PRIMARY KEY (`stateID`)
) ENGINE=MyISAM AUTO_INCREMENT=273 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `transfer`
--

DROP TABLE IF EXISTS `transfer`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `transfer` (
  `transferID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seqProjectID` int(10) unsigned DEFAULT NULL,
  `transLocation` varchar(500) DEFAULT NULL,
  `sftpAccountID` int(10) unsigned DEFAULT NULL,
  `dataStatus` enum('origin','sftp') DEFAULT NULL,
  PRIMARY KEY (`transferID`)
) ENGINE=MyISAM AUTO_INCREMENT=199 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `typeLinker`
--

DROP TABLE IF EXISTS `typeLinker`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `typeLinker` (
  `linkID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` enum('masterProject','sequencing','flow cytometry','analysis') DEFAULT NULL,
  `parentID` int(10) unsigned DEFAULT NULL,
  `childID` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`linkID`)
) ENGINE=MyISAM AUTO_INCREMENT=206 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2015-09-25 12:27:47
