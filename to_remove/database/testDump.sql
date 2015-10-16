-- MySQL dump 10.13  Distrib 5.5.32, for Linux (x86_64)
--
-- Host: localhost    Database: coreInSys
-- ------------------------------------------------------
-- Server version	5.5.32

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
) ENGINE=InnoDB AUTO_INCREMENT=13 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `IDtagLibs`
--

LOCK TABLES `IDtagLibs` WRITE;
/*!40000 ALTER TABLE `IDtagLibs` DISABLE KEYS */;
INSERT INTO `IDtagLibs` VALUES (1,'Illumina','1','ATCACG'),(2,'Illumina','2','CGATGT'),(3,'Illumina','3','TTAGGC'),(4,'Illumina','4','TGACCA'),(5,'Illumina','5','ACAGTG'),(6,'Illumina','6','GCCAAT'),(7,'Illumina','7','CAGATC'),(8,'Illumina','8','ACTTGA'),(9,'Illumina','9','GATCAG'),(10,'Illumina','10','TAGCTT'),(11,'Illumina','11','GGCTAC'),(12,'Illumina','12','CTTGTA');
/*!40000 ALTER TABLE `IDtagLibs` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `customer`
--

DROP TABLE IF EXISTS `customer`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `customer` (
  `customerID` int(11) NOT NULL DEFAULT '0',
  `name` varchar(40) DEFAULT NULL,
  `email` varchar(50) DEFAULT NULL,
  `tel` varchar(30) DEFAULT NULL,
  PRIMARY KEY (`customerID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `customer`
--

LOCK TABLES `customer` WRITE;
/*!40000 ALTER TABLE `customer` DISABLE KEYS */;
/*!40000 ALTER TABLE `customer` ENABLE KEYS */;
UNLOCK TABLES;

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
) ENGINE=InnoDB AUTO_INCREMENT=7 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `laneData`
--

LOCK TABLES `laneData` WRITE;
/*!40000 ALTER TABLE `laneData` DISABLE KEYS */;
INSERT INTO `laneData` VALUES (1,1,5,0.5,0,0.5,'NULL',0),(2,2,5,0.5,0,0.5,'NULL',0),(3,3,6,0.5,0,0.5,'NULL',0),(4,1,7,0.5,0,0.5,'NULL',0),(5,2,7,0.5,0,0.5,'NULL',0),(6,3,8,0.5,0,0.5,'NULL',0);
/*!40000 ALTER TABLE `laneData` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `masterProject`
--

DROP TABLE IF EXISTS `masterProject`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `masterProject` (
  `masterProjectID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `projectType` enum('sequencing','flow cytometry','other') DEFAULT NULL,
  `status` varchar(20) NOT NULL,
  `projectName` varchar(20) DEFAULT NULL,
  `localProjectID` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`masterProjectID`)
) ENGINE=InnoDB AUTO_INCREMENT=5 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `masterProject`
--

LOCK TABLES `masterProject` WRITE;
/*!40000 ALTER TABLE `masterProject` DISABLE KEYS */;
INSERT INTO `masterProject` VALUES (1,'sequencing','analysis','Piranha',8),(2,'sequencing','analysis','Hammerhead',6),(3,'sequencing','analysis','Conger eel',5),(4,'sequencing','analysis','Stingray',7);
/*!40000 ALTER TABLE `masterProject` ENABLE KEYS */;
UNLOCK TABLES;

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
  `analysisID` int(10) unsigned DEFAULT NULL,
  `adaptorSequence` varchar(200) DEFAULT NULL,
  PRIMARY KEY (`sampleID`)
) ENGINE=InnoDB AUTO_INCREMENT=31 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `sampleData`
--

LOCK TABLES `sampleData` WRITE;
/*!40000 ALTER TABLE `sampleData` DISABLE KEYS */;
INSERT INTO `sampleData` VALUES (1,'A1',1,1,'AAAAAA',1,NULL),(2,'A2',2,1,'GAAAAA',1,NULL),(3,'A3',3,1,'AGAAAA',1,NULL),(4,'A4',4,1,'AAGAAA',1,NULL),(5,'A5',5,1,'AAAGAA',1,NULL),(6,'A1',1,2,'AAAAAA',1,NULL),(7,'A2',2,2,'GAAAAA',1,NULL),(8,'A3',3,2,'AGAAAA',1,NULL),(9,'A4',4,2,'AAGAAA',1,NULL),(10,'A5',5,2,'AAAGAA',1,NULL),(11,'A1',1,3,'AAAAAA',1,NULL),(12,'A2',2,3,'GAAAAA',1,NULL),(13,'A3',3,3,'AGAAAA',1,NULL),(14,'A4',4,3,'AAGAAA',1,NULL),(15,'A5',5,3,'AAAGAA',1,NULL),(16,'A1',1,4,'AAAAAA',1,NULL),(17,'A2',2,4,'GAAAAA',1,NULL),(18,'A3',3,4,'AGAAAA',1,NULL),(19,'A4',4,4,'AAGAAA',1,NULL),(20,'A5',5,4,'AAAGAA',1,NULL),(21,'A1',1,5,'AAAAAA',1,NULL),(22,'A2',2,5,'GAAAAA',1,NULL),(23,'A3',3,5,'AGAAAA',1,NULL),(24,'A4',4,5,'AAGAAA',1,NULL),(25,'A5',5,5,'AAAGAA',1,NULL),(26,'A1',1,6,'AAAAAA',1,NULL),(27,'A2',2,6,'GAAAAA',1,NULL),(28,'A3',3,6,'AGAAAA',1,NULL),(29,'A4',4,6,'AAGAAA',1,NULL),(30,'A5',5,6,'AAAGAA',1,NULL);
/*!40000 ALTER TABLE `sampleData` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `seqRun`
--

DROP TABLE IF EXISTS `seqRun`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `seqRun` (
  `seqRunID` varchar(10) NOT NULL DEFAULT '',
  `flowcellID` varchar(20) DEFAULT NULL,
  `startDate` date DEFAULT NULL,
  `completionDate` date DEFAULT NULL,
  `genomicsLead` varchar(20) DEFAULT NULL,
  `dataLocation` varchar(40) DEFAULT NULL,
  `indexTagCycles` tinyint(4) DEFAULT NULL,
  `readCycles` tinyint(4) DEFAULT NULL,
  PRIMARY KEY (`seqRunID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `seqRun`
--

LOCK TABLES `seqRun` WRITE;
/*!40000 ALTER TABLE `seqRun` DISABLE KEYS */;
INSERT INTO `seqRun` VALUES ('P211',NULL,NULL,NULL,NULL,NULL,NULL,NULL),('P222',NULL,NULL,NULL,NULL,NULL,NULL,NULL);
/*!40000 ALTER TABLE `seqRun` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `seqProject`
--

DROP TABLE IF EXISTS `seqProject`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `seqProject` (
  `seqProjectID` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seqProjectName` varchar(20) DEFAULT NULL,
  `masterProjectID` int(10) unsigned DEFAULT NULL,
  `seqRunID` varchar(10) DEFAULT NULL,
  `customerID` int(11) DEFAULT NULL,
  PRIMARY KEY (`seqProjectID`)
) ENGINE=InnoDB AUTO_INCREMENT=9 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `seqProject`
--

LOCK TABLES `seqProject` WRITE;
/*!40000 ALTER TABLE `seqProject` DISABLE KEYS */;
INSERT INTO `seqProject` VALUES (5,'AlProj',3,'P211',0),(6,'AlPrxxoj',2,'P211',0),(7,'AlProjb',4,'P222',0),(8,'AlP11xxoj',1,'P222',0);
/*!40000 ALTER TABLE `seqProject` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2013-07-11 16:17:37
