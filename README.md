# Biocompute-DM

## Overview
BioCompute-DM is a software based infrastructure for data management and processing. It facilitates project-centric tracking through centralised document and data storage, as well as provides an intuitive ‘computational process’ execution interface.
BioCompute-DM utilises a small but powerful work-flow model to cover common data management and analysis use-cases with a core focus towards next generation sequencing problems. 
BioCompute-DM work-flows are augmented through plugin-able pipelines, where each pipeline provides unique intelligence for specific data processing problems. BioCompute-DM pipelines require a small structure be wrapped around traditional bioinformatic solutions which, in turn, imparts comprehensive options systems and automated tracking solutions to each pipeline providing end-users with rapid, robust and flexible data analysis.
BioCompute-DM is built upon a database infrastructure that is agnostic to the data types and processing methods. It provides a fully auditable layer for building analytics intelligence whilst allowing the user to fully manage their projects. 
Built in to the infrastructure is the concept of users and clients, which allows data sharing between groups whilst imparting two permission levels (those with and without permission to execute new analyses).
	
## Source Code & Implementation

### Core Frameworks:

#### Flask Web page:
BioCompute-DM uses the flask microframework to serve it’s web pages. Flask is based on python and provides a series of pre-built modules that promote rapid site development. The following modules are used within BioCompute-DM:

* Flask-Script
* PyMySQL
* flask-sqlalchemy
* alembic
* flask-migrate
* flask-login
* flask-mail
* flask-wtf
* flask-bootstrap
* jsonschema

The structure of the flask applications source code helps delimit core concepts for the application. Predominantly we have several subdirectories which represent thematic sections of the applications content. These are usually mostly self contained; they have a sub directories containing the html page templates and an optional helper sub directory that may contain logic items. Additionally the thematic sections traditionally have 3 core files termed: forms, models and views. Forms contain any web page forms that are used within their section, this file is optional to the themed sections. Views contains the access point for each web page end point within the website (relevant to that section) and models contains any corresponding database items. These three concepts are normally tightly linked.
Also within the source directory is a subdirectory named static which contains all of the additional files that must be served from the web address (and is also used as a fake object server for allowing items to be downloaded). This directory includes the embedded java web applet for file transfers.
Finally the templates directory contains a series of base html items that are usually extended upon by the more specific subdirectories discussed above (i.e. they may provide common views etc).

#### MySQL Database:
The underlying database used by BioCompute-DM is a MySQL database that is mapped within the flask application (using the pymysql flask-sqlalchemy modules). A schematic diagram of the database can  be found below. The model contains three core information layers: the pipeline layer, the data layer and the management layer.
The flask modules also include flask-migrate. This plugin allows migration scripts to be generated per database iteration to ensure that live data is never damaged. Database upgrade scripts can be created once the flask models have been updated. These migration scripts can then be applied to a live database.

### Configurable Frameworks:
BioCompute-DM also contains a high degree of modularity and configurability to maximise it’s flexibility when being deployed. To achieve this, several of the core frameworks have been abstracted sufficiently to allow interoperability between systems. Generally this flexibility is managed through a config.py file that must be created (according to the config_template.py file). Below is a list of their values as well as a description of their function.

* SECRET_KEY: Key used to encrypt sessions and cookies
* ERROR_LOG_LOCATION: The location to store the error log, it must be accessible by the biocompute-dm user
* MAIL_USERNAME: The username for your mail account
* MAIL_PASSWORD: The password for your mail account
* MAIL_DEFAULT_SENDER: The sender address for your mail
* MAIL_SERVER: The mail server for your mail account
* MAIL_PORT: The accessible port for your mail account
* MAIL_USE_TLS: (Refer to flask-mail)
* MAIL_USE_SSL: (Refer to flask-mail)
* DATABASE_USERNAME: The username for your database account
* DATABASE_PASSWORD: The password for your database account
* DATABASE_LOCATION: The network address of your database (suggested localhost)
* DATABASE_NAME: The name of your database
* SITE_ADMIN_USERNAME: The username of the site admin
* SITE_ADMIN_PASSWORD: The password of the site admin
* SITE_ADMIN_EMAIL: The email of the site admin
* SITE_GROUP_PASSWORD: The password for the site admin group
* HPC_DEBUG: Flag to use alternative job submission script (fake_submit_job.sh)
* LOCAL_WEBSERVER_PORT: The port your webserver runs on
* NETWORK_PATH_TO_WEBSERVER_FROM_HPC: The ip address to connect to hpc resource (which will be ssh’d to by the application)
* NETWORK_PATH_TO_HPC_FROM_WEBSERVER: The ip address to reach the webserver from the hpc resource
* HPC_USERNAME: The hpc username (should have passwordless permissions such as ssh keys)
* HPC_ROOT_PATH: The path to the data on the hpc
* WEBSERVER_ROOT_PATH: The path to the data on the webserver
* TEMP_DIRECTORY: Must be accessible to the biocompute-dm user
* SFTP_USER_ROOT_PATH: Directory containing the sftp root path (where all user sftp directories will be placed and chrooted) – must be owned by root
* MANAGEMENT_SCRIPTS_PATH_AFTER_RELATIVE_ROOT: The path following either the hpc root path or the webserver root path to reach management scripts
* PIPELINE_DATA_PATH_AFTER_RELATIVE_ROOT: As above but for pipeline data
* SUBMISSION_DATA_PATH_AFTER_RELATIVE_ROOT: As above but for submission data
* SAMPLE_DATA_PATH_AFTER_RELATIVE_ROOT: As above but for sample data
* PROJECT_DATA_PATH_AFTER_RELATIVE_ROOT: As above but for project data
* REFERENCE_DATA_PATH_AFTER_RELATIVE_ROOT: As above but for reference data

#### Underlying Operating System:
Most of the core functionality of the system (with regards to underlying operating system operations) can be amended by replacing or modifying the scripts found under biocompute-dm/scripts. These usually include the admin, cleanup, and io subdirectories.

#### Compute Execution:
In the same vein as above the underlying hpc infrastructure can also be changed by modifying the files in the jobs subdirectory to better reflect specific installations. 

#### Data Focus:
Data types are very well abstracted in BioCompute-DM but the core functionality of any data processing / analysis can be modified by creating / modifying the pipelines present in the pipeline subdirectory. 

### Additional Notes
One of the most important files to investigate in order to understand how some of the webserver/hpc logic breaks up is the utils.py function get_path as this is called often to find side specific urls.

## Installation and Routine Tasks

### Installation
With a LAMP(y) stack base (Python v3), create a group called sftpusers, and a user called biocompute-dm (add the biocompute-dm user to sftp users).
Ensure your biocompute-dm user has ssh access without password to your HPC environment. 
Ensure the following packages are installed:

* 7z
* unrar
* unzip

Clone the biocompute-dm bitbucket repo into /var/www (or equivalent)
Create a python3 virtual environment within the newly created /var/www/biocompute-dm

* python3 -m venv flask_environment

Install the following packages using pip:

* flask_environment/bin/pip install flask
* flask_environment/bin/pip install Flask-Script
* flask_environment/bin/pip install PyMySQL
* flask_environment/bin/pip install flask-sqlalchemy
* flask_environment/bin/pip install alembic
* flask_environment/bin/pip install flask-migrate
* flask_environment/bin/pip install flask-login
* flask_environment/bin/pip install flask-mail
* flask_environment/bin/pip install flask-wtf
* flask_environment/bin/pip install flask-bootstrap
* flask_environment/bin/pip install jsonschema

Create / modify an apache2 conf (in sites-enabled) with the following properties (adjustments to your needs may be necessary). The core consideration here is that the apache execution user is biocompute-dm for this website.
Ensure your HPC parent environment (e.g. head nodes) have access to your apache2 port 80 (even if you only wish to use https, internally port 80 must still be available).
Setup the following properties to enable sftp to your server (In /etc/ssh/sshd_config):

* /#Subsystem sftp /usr/libexec/openssh/sftp-server
* Sybsystem sftp internal-sftp
* 
* /#Add the following to the end of the file:
* Match Group sftpusers
*   ChrootDirectory %h
*   ForceCommand internal-sftp
*   AllowTcpForwarding no

Ensure that all of the scripts in your biocompute-dm install have sudo rites and ensure that all of these files are owned by route and have 700 permissions:

* biocompute-dm ALL=(ALL) NOPASSWD: [path to script]

Create a custom config.py according to your needs. Ensure that all of the relevant files are accessible on both the HPC and the webserver environments (mapped drives etc). You may need to copy your script files into a shared location.
Create a database according to the properties declared in your config.py with correct permissions for localhost access by the biocompute-dm user.
Execute the first time setup of your database using the following command (relative to the biocompute-dm folder – as biocompute-dm):

* flask_environment/bin/python manage.py update

### Routine Tasks

#### Management
Several key database management scripts are provided to help manage live servers (each executed as the command above):

* clean FORCE:
	wipe the db to default, restore only the admin user

* migrate:
    used to generate migration scripts for production from dev environmets

* update:
    use the available migration scripts to update production servers

#### Add New Pipelines
New pipelines are added by placing their folders (see pipeline development documentation for details) in your scripts/pipelines directory. From there, access the site administration admin panel and choose refresh pipelines.

#### Add New Reference File
This is performed as above but the folders are replaced within your reference directory. Instead of refreshing pipelines, use the refresh reference files button in the administration panel, Reference files should use the following json template:

    {
       “name”: “real name of file”,
       “description”: “displayable description”,
       “version”
    }

## Pipeline Development
Pipelines are integrated into the BioCompute-DM framework through a common template entry point. The underlying assumptions for pipeline integration are as follows:

* Pipeline contents (scripts etc) are nested in one common folder.
* The folder name is the same as the filled template (within the folder)

Pipelines are made up of modules (minimum of 1) that represents core steps with the pipeline work flow. Each module is represented by a bash executable (which in turn can call any language or software available in the HPC environment). Each module is provided with a core set of environmental variables that allow it to appropriately understand its current location as well as the location of important resources.

The following environmental variables are provided as standard:

* USERNAME: The current username within the HPC environment
* HPC_IP: The IP address to ssh to back to parent HPC infrastructure
* TICKET: A unique ticket id for your current execution
* PIPELINE_SOURCE_DIRECTORY: The directory of the pipeline executors, allows resources to be accessed
* PIPELINE_OUTPUT_DIRECTORY: The directory for pipeline specific (i.e. non sample) outputs
* MODULE_OUTPUT_DIRECTORY: The directory for module specific outputs – traditionally used for logs etc
* WORKING_DIRECTORY: The directory in which you should be executing logic
* SERVER: The IP address to contact the webserver – usually only available when ssh’d back to the HPC_IP (see helpful code example in reference).

Modules are also provided a link to a small csv (corresponding to the environmental variable SAMPLE_CSV) file that contains  information about the source data locations. This csv file is structured as follows:

* Data Item Name, Data Item Location, Data Item Output, Empty

The template also allows the modules to request custom environmental variables, denote whether or not these variables are necessary or optional as well as provide default implementations. These custom variables are further described in the template explanation.

#### Pipeline Template:
      {
        “name”: “pipeline name”,
        “description”: “pipeline description”,
        “pipeline type”: “see comment 1”,
        “author”: “author, author2”,
        “version”: “semver version”,
        “regex_type”: “see comment 2”,
        “file_regex”: [“see comment 3”],
        “documentation_file_name”: “see comment 4”.
        “modules”: [
          { // see comment 5
            “name”: “module name”,
            “description”: “module description”,
            “executor”: “see comment 6”,
            “options”: [
              { // see comment 7
                “display_name”: “Human readable option text”,
                “parameter_name”: “variable name to inject to executor”,
                “default_value”: “value to use if defaults are chosen”,
                “user_interaction_type”: “see comment 8”
                “necessary”: see comment 9,
                “description”: “Descriptive information for users”
              }
            ]
          }
        ]
      }

##### Comments:

1. Pipelines can be of the following types (string): I, II, III. This allows granularity over when the pipelines can be executed. I type pipelines (or pre-processors) are responsible for loading raw data into BioCompute-DM and producing samples as an end point. These samples represent the base unit of operation across all subsequent processing. II and III (processing and post-processing) are currently available interchangeably but in future could be used to further enforce logical process ordering.
2. Pipelines use regular expressions to check that the data contains certain file types (to prevent pipelines being initiated on data that cannot be processed). The regex_type flag corresponds to one of two variables: “OR” or “AND”. The former indicates that any of the regexes provided in the file_regex array be positive whilst the latter indicates that all must be present.
3. Here you may provide a comma separated list (with quotes) of java compliant (i.e. double escape - \\ - for literals) regexes to match against your source data.
4. A file name for a file that is present within your pipeline directory that links to your documentation. This will be made available from the web page.
5. This is an object array of modules that will be processed in order when the pipeline is executed. Multiple modules are added by adding a comma after the first module’s curly braces and then opening another set.
6. This is the name of the shell script that you wish to execute. The shell script must be present within the same folder as the json.
7. This is an object array of options that can be provided to the corresponding module. Multiple options are added in the same way as multiple modules.
8. User interaction types may take 5 forms: “boolean”, “string”, “reference”, “file”, and “enum”. Each of these corresponds to the type of information expected from the user. Below is a list of what each means:
    * “boolean”: provides the user with a tick box, which in turn provides the module execution environment with a “True” or “False” string
    * “string”: provides the user with a text box in which to type information. This is provided to the module environment with one exception, all commas are replaced with the following flag: %%__%% - this must be handled by the module executors accordingly.
    * “reference”: Provides the user with a drop down list of files available in the reference database of BioCompute-DM. 
    * “enum”: Provides the user with a set of text based options. The options are generated from the value placed in the “default_value” field where commas are separated into different options and the first value listed represents the actual default.
    * “file”: Allows the user to upload a file which will be provided to the module.
9. The necessary flag, when sets to true, ensures that the user must provide a value for this option (i.e. ignore any default value)

#### Reference:
The following excerpt is a useful pattern for handling when a pipeline has errored:

* ssh ${USERNAME}@${HPC_IP} << EOF
* curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
* EOF

For reference, the pipeline template is validated against the following json schema rules:

    {
      "type": "object",
      "properties": {
        "name": {
          "type": "string"
        },
        "description": {
          "type": "string"
        },
        "pipeline_type": {
          "enum": ["I", "II", "III"]
        },
        "author": {
          "type": "string"
        },
        "version": {
          "type": "string"
        },
        "regex_type": {
          "enum": ["AND", "OR"]
        },
        "file_regex": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "minItems": 1,
          "uniqueItems": true
        },
        "documentation_file_name": {
          "type": "string"
        },
        "modules": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "name": {
                "type": "string"
              },
              "description": {
                "type": "string"
              },
              "executor": {
                "type": "string"
              },
              "options": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "display_name": {
                      "type": "string"
                    },
                    "parameter_name": {
                      "type": "string"
                    },
                    "default_value": {
                      "type": "string"
                    },
                    "user_interaction_type": {
                      "enum": ["boolean", "string", "reference", "file", "enum"]
                    },
                    "necessary": {
                      "type": "boolean"
                    },
                    "description": {
                      "type": "string"
                    }
                  },
                  "additionalProperties": false,
                  "required": [
                    "display_name",
                    "parameter_name",
                    "default_value",
                    "user_interaction_type",
                    "necessary",
                    "description"
                  ]
                }
              }
            },
            "additionalProperties": false,
            "required": [
              "name",
              "description",
              "executor",
              "options"
            ]
          }
        }
      },
      "additionalProperties": false,
      "required": [
        "name",
        "description",
        "pipeline_type",
        "author",
        "version",
        "regex_type",
        "file_regex",
        "documentation_file_name",
        "modules"
      ]
    }

