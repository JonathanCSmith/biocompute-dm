# Biocompute-DM
## Repository for the Biocompute Data Manager

#### Note, this branch is still very WIP. Viability of the software on a commit by commit basis cannot be guaranteed!

### Ubuntu Installation (mileage will vary for other OSs)
- Setup a machine with a Lamp(y) stack
    - Linux, apache2, mysql, python3

- Create a group called sftpusers, a user called biocompute-dm and add biocompute-dm to the group sftpusers

- Clone branch containing the flask builds onto the machine to /var/www (or equivalent)

- Within the repository folder create a virtual environment:

        - python3 -m venv flask_environment

- OR (if this fails):

        - python3 -m venv --without-pip flask_environment
        - source flask_environment/bin/activate
        - curl https://bootstrap.pypa.io/get-pip.py | python
        - deactivate

- Then install flask packages using the following commands:

        - flask_environment/bin/pip install flask
        - flask_environment/bin/pip install Flask-Script
        - flask_environment/bin/pip install PyMySQL
        - flask_environment/bin/pip install flask-sqlalchemy
        - flask_environment/bin/pip install alembic
        - flask_environment/bin/pip install flask-migrate
        - flask_environment/bin/pip install flask-login
        - flask_environment/bin/pip install flask-mail
        - flask_environment/bin/pip install flask-wtf
        - flask_environment/bin/pip install flask-bootstrap
        - flask_environment/bin/pip install jsonschema
        
- Modify your apache2 installation according to the best practices listed below
  * It is always assumed that script execution on the webserver is performed using the user 'biocompute-dm'
  * This administrative user must be a member of the group sftpusers

- SFTP Setup
1. Openssh sftp is allowed and configured correctly using the /etc/ssh/sshd_config file
  * In the middle of the file comment out the first line and add the second:
  
        #Subsystem       sftp    /usr/libexec/openssh/sftp-server
        Sybsystem sftp internal-sftp
        
  * At the end of the file add the following information:
    
        Match Group sftpusers
           ChrootDirectory %h
           ForceCommand internal-sftp
           AllowTcpForwarding no
           
2. The sftpusers group exists in the system
3. The sftp bash scripts have the correct permissions to be executed as sudo
  * Type sudo visudo
  * Below the line `%sudo ALL=(ALL:ALL) ALL` add the followings for each sftp (in your sftp folder) bash script:
        
        biocompute-dm ALL=(ALL) NOPASSWD: [path to script]
        
  * Ensure that the files are owned by root and have 700 permissions! (this prevents nefarious execution)

- Create your own config.py by copying the template and inserting the relevant information.

- Follow database instructions for "New database schema"
  * Note: all commands within manage.py should be called as your biocompute-dm user (or the owner of the error.log directory)
  
- Install the following unpacking programs:
  * 7z
  * unrar
  * unzip

### Database Commands - TODO: This information is outdated - update!!
#### New database schema

- Create a mysql db according to your config.py
- Setup a user with properties that corresponds to your config.py for accessing the db above (local only)
- as biocompute-dm user in the website root: flask_environment/bin/python manage.py update

### Pipeline Patterns - Rolling your own

TODO Overview text here

#### Pipeline Types - TODO: This information is outdated - update!!
There are three pipeline types that can be created in Biocompute-DM, respectfully represented by their roman numerals. 

* Type I plugins deconvolute information into separate samples. This pipeline type is markedly different from its counterparts
  as it directly informs the database of additional information by analysing the path on which the pipeline's output
  was stored. I.e. within .../pipeline_instance_id/... the pipeline informs that database about separate samples by generatring folders
  for them based on internally assigned sample ids and each folder must contain a sample.json file (see below). There will be no namespace conflict other than that used by the pipeline as 
  these names are later translated into internal GUIDs.
  
        The template sample.json that must be filled in for each sample folder (nothing is enforced, it just helps the users)
        {
          "name": "sample_name",
          "description": "sample_description",
          "notes": "sample_notes_500_characters",
          "internal_reference_id": "internal_ref",
          "client_reference_id": "client_ref"
        } TODO: We don't have this anymore - remove!!!
  
* Type II plugins work on a group of samples produced by the same Type I pipeline instance. This normally represents the
 processing stage.
 

* Type III plugins work on groups of sample groups, normally attributed to a project and may contain cross-sample, cross-batch
 type controls.
 
#### Pipeline Properties


#### Pipeline module options


#### Pipeline Executor Environmental Variables

 
#### Example Pipeline - TODO: This information is outdated - update!!
Below is an example pipeline. Things to note:

* pipeline_type is an enum of I, II or III
* the modules are executed in the order that they are listed
* parameter name represents the flag that will be provided to your executor with the value selected by the user (in the form parameter=value)
* user_interaction_type is an enum of string, boolean, library


    {
      "name": "example_pipeline",
      "description": "this is a test - note it must reside in a folder with the same name! i.e. example_pipeline",
      "pipeline_type": "I",
      "author": "me",
      "version": "0.0.0",
      "modules": [
        {
          "name": "module_1",
          "description": "step 1",
          "executor": "somefile",
          "options": [
            {
              "display_name": "Cool Selection!",
              "parameter_name": "cs",
              "default_value": "one",
              "user_interaction_type": "string"
            }
          ]
        },
        {
          "name": "module_2",
          "description": "step 2, note the executor for this module must be present in the same location as this file",
          "executor": "somefile2",
          "options": [
            {
              "display_name": "Cool Selection!",
              "parameter_name": "csb",
              "default_value": "",
              "user_interaction_type": "boolean"
            },
            {
              "display_name": "Library Selection!",
              "parameter_name": "l",
              "default_value": "",
              "user_interaction_type": "library"
            }
          ]
        }
      ]
    }
    
### Reference Files

