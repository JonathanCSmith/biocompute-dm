# BioCIS
## Repository for the BioCIS pipeline

#### Note, this branch is still very WIP. Viability of the software on a commit by commit basis cannot be guaranteed!

### Ubuntu Installation (mileage will vary for other OSs)
Setup a machine with a Lamp(y) stack

Clone branch containing the flask builds onto the machine to /var/www (or equivalent)

Within the biocis folder create a virtual environment:

- python3 -m venv flask

OR (if this fails):

- python3 -m venv --without-pip flask
- source flask/bin/activate
- curl https://bootstrap.pypa.io/get-pip.py | python
- deactivate

Then install flask packages using the following commands:

- flask/bin/pip install flask
- flask/bin/pip install flask-sqlalchemy
- flask/bin/pip install sqlalchemy-migrate
- flask/bin/pip install flask-login
- flask/bin/pip install flask-mail
- flask/bin/pip install flask-whooshalchemy
- flask/bin/pip install flask-wtf
- flask/bin/pip install flask-babel
- flask/bin/pip install guess_language
- flask/bin/pip install flipflop
- flask/bin/pip install coverage
- flask/bin/pip install PyMySQL
- flask/bin/pip install flask-excel
- flask/bin/pip install flask-pyexcel-xls
- flask/bin/pip install flask-pyexcel-xlsx
- flask/bin/pip install flask-bootstrap
- flask/bin/pip install jsonschema

Modify your apache2 installation according to the best practices listed below

Follow database instructions for a "first time setup"

Generate the workspace directories by creating a folder called link that contains symlinks to the following:
- 'Investigations' <- Where you wish to house your investigation associated information
- 'Pipelines' <- Where existing / new pipelines are placed for execution. Note pipeline patterns are documented below.
- 'Samples' <- Where sample data will sit after it has been pre-processed by pipeline type is
- 'Raw' <- Where the raw data will sit after submission to biocis (either through automated or manual means)

Finally, create your own config.py by copying the template and inserting the relevant information.
This information will need to contain the name of your database as well as user information with rights to access said database.

### Database Commands
#### New database schema

- Create a mysql db
- run db_create.py
- run db_migrate.py

#### First time setup

- run db_create.py
- run db_upgrade.py

#### Schema update

- run db_migrate.py

### Apache2 Best Practices

- TODO: Bug me if you need this asap


### Pipeline Patterns - Rolling your own

#### Pipeline Types
There are three pipeline types that can be created in BioCIS, respectfully represented by their roman numerals. 

Type I plugins deconvolute information into separate samples. This pipeline type is markedly different from its counterparts
  as it directly informs the database of additional information by recursively analysing the path on which the pipeline's output
  was stored. I.e. within .../pipeline_instance_id/... the pipeline informs that database about separate samples by generatring folders
  for them based on internally assigned sample ids. There will be no namespace conflict other than that used by the pipeline as 
  these names are later translated into internal GUIDs.
  
Type II plugins works on a group of samples produced by the same Type I pipeline instance. This normally represents the
 processing stage.
 
Type III plugins work on groups of sample groups, normally attributed to a project and may contain cross-sample, cross-batch
 type controls.

- TODO

