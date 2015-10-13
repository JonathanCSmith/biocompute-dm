# biocis
## Repository for the BioCIS pipeline

### Installation
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

Finally follow database instructions for a "first time setup"

### Database Commands

#### New database schema

- run db_create.py
- run db_migrate.py

#### First time setup

- run db_create.py
- run db_upgrade.py

#### Schema update

- run db_migrate.py



### To Setup (DEPRECATED):

  Setup a machine with a Lamp(y) stack
  
  Clone branch onto machine for deployment
  
  Move all files (except db schema) into apache2 website location
  
  Build using edit.py in the build folder
  
  Link apache2 website root directory to /biocis/
  
  Link apache2 directory root (for the website) to /biocis/html/index.html
  
  Deploy mysql schema into db
  
  Change configuration files to reflect db propertes
  
  Profit?
