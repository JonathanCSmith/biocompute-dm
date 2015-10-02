__author__ = 'jon'

from app import app, db, forms, login_required
from flask import render_template, request, flash, redirect, url_for, g, session
from flask.ext.login import login_user, logout_user, current_user


@app.route("/")
@app.route("/index")
def index():
    return render_template("index.html", title="Home", user=g.user)


@app.route("/register", methods=["GET", "POST"])
def register():
    if g.user is not None and g.user.is_authenticated:
        return redirect(url_for("index"))

    form = forms.RegisterForm()
    if request.method == "GET":
        return render_template("register.html", form=form)

    else:
        if not form.validate():
            return render_template("register.html", form=form)

        else:
            from app.models import User
            user = User(str(form.username.data), str(form.password.data), str(form.email.data.lower()))
            db.session.add(user)
            db.session.commit()
            flash("User successfully registered!")
            # session["remember_me"] = form.remember_me.data
            # app.login_user(user, remember=form.remember_me.data)
            login_user(user)

        return redirect(url_for("index"))


@app.route("/login", methods=["GET", "POST"])
def login():
    if g.user is not None and g.user.is_authenticated:
        return redirect(url_for("index"))

    form = forms.LoginForm()
    if request.method == "GET":
        return render_template("login.html", form=form)

    else:
        if not form.validate_on_submit():
            return render_template("login.html", form=form)

        else:
            from app.models import User
            user = User.query.filter_by(username=str(form.username.data)).first()
            # session["remember_me"] = form.remember_me.data
            # app.login_user(user, remember=form.remember_me.data)
            flash("Successfully logged in!")
            login_user(user)
            return redirect(url_for("index"))


@app.route("/logout")
@login_required("ANY")
def logout():
    logout_user()
    return redirect(url_for("index"))


@app.route("/projects", methods=["GET", "POST"])
@app.route("/projects/<int:page>", methods=["GET", "POST"])
@login_required()
def projects(page=1):
    from app.models import MasterProject
    projects = MasterProject.query.paginate(page=page, per_page=20)
    flash("Projects page is still in development")
    return render_template("projects.html", page=page, projects=projects)


@app.route("/project/<string:name>|<int:id>", methods=["GET", "POST"])
@login_required()
def project(name="", id=-1):
    if name == "" or id == -1:
        flash("Incorrect arguments for query provided.")
        return redirect(url_for("index"))

    from app.models import MasterProject
    project = MasterProject.query.filter_by(masterProjectID=id, projectName=name).first()

    return redirect(url_for("index"))


@app.route("/sequencing_runs")
@login_required("ANY")
def sequencing_runs():
    flash("Sequencing runs page is still in development")
    return redirect(url_for("index"))


@app.route("/input")
@login_required("ANY")
def input():
    flash("Input page is still in development")
    return render_template("input.html")


@app.route("/new_project", methods=["GET", "POST"])
@login_required("ANY")
def new_project():
    flash("New project page is still in development")
    form = forms.NewProjectForm()
    if request.method == "GET":
        return render_template("new_project.html", form=form)

    else:
        if form.validate_on_submit():
            from app.models import MasterProject
            p = MasterProject(str(form.project_name.data), str(form.project_lead.data))
            db.session.add(p)
            db.session.commit()
            flash("Project successfully registered!")
            # session["remember_me"] = form.remember_me.data
            # app.login_user(user, remember=form.remember_me.data)

        return render_template("new_project.html", form=form)


@app.route("/input_sequencing_run")
@login_required()
def input_sequencing_run():
    flash("New sequencing run page is still in development")
    return redirect(url_for("index"))


@app.route("/input_sequencing_project")
@login_required()
def input_sequencing_project():
    flash("New sequencing project page is still in development")
    return redirect(url_for("index"))

# @login_required to secure
