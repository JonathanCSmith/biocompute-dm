import os
import flask.ext.excel as excel
import pyexcel.ext.xls
import pyexcel.ext.xlsx
import pyexcel.ext.ods

from app import app, db, forms, login_required
from flask import render_template, request, flash, redirect, url_for, g, jsonify
from flask.ext.login import login_user, logout_user

__author__ = 'jon'


@app.route("/")
@app.route("/index")
def index():
    return render_template("index.html", title="Home", user=g.user)


@app.route("/empty")
def empty():
    return render_template("empty.html")


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
            flash("User successfully registered!", "success")
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
            flash("Successfully logged in!", "success")
            login_user(user)
            return redirect(url_for("index"))


@app.route("/logout")
@login_required("ANY")
def logout():
    logout_user()
    return redirect(url_for("index"))


@app.route("/investigations", methods=["GET", "POST"])
@app.route("/investigations/<int:page>", methods=["GET", "POST"])
@login_required("ANY")
def investigations(page=1):
    from app.models import Investigation
    i = Investigation.query.paginate(page=page, per_page=20)
    return render_template("investigations.html", page=page, investigations=i)


@app.route("/new_investigation", methods=["GET", "POST"])
@login_required("ANY")
def new_investigation():
    form = forms.NewInvestigationForm()
    if request.method == "GET":
        return render_template("new_investigation.html", form=form)

    else:
        if form.validate_on_submit():
            from app.models import Investigation
            p = Investigation(str(form.investigation_name.data), str(form.investigation_lead.data))
            db.session.add(p)
            db.session.commit()
            flash("Investigation successfully registered!", "info")
            return redirect(url_for("investigations"))

        return render_template("new_investigation.html", form=form)


@app.route("/investigation/<string:name>|<int:iid>", methods=["GET", "POST"])
@login_required("ANY")
def investigation(name="", iid=-1):
    if name == "" or iid == -1:
        flash("Incorrect arguments for query provided.", "error")
        return redirect(url_for("index"))

    from app.models import Investigation, Document
    i = Investigation.query.filter_by(investigation_id=iid, investigation_name=name).first()

    return render_template("investigation.html", investigation=i, projects=[])


@app.route("/add_document/<int:iid>", methods=["GET", "POST"])
@login_required("ANY")
def add_document(iid=-1):
    if iid == -1:
        flash("Incorrect arguments for query provided", "error")
        return redirect(url_for("index"))

    form = forms.AddDocumentForm()
    if request.method == "POST":
        if form.validate_on_submit():
            from app.models import Investigation, Document
            i = Investigation.query.filter_by(investigation_id=iid).first()
            if not i:
                flash("Attempted to add a document to an investigation without an investigation parent.", "error")
                return redirect(url_for("investigations"))

            directory = i.validate_investigation_directory()
            if directory is None:
                flash("Could not create investigation directory.", "error")
                return redirect(url_for("investigations"))

            from werkzeug.utils import secure_filename
            filename = secure_filename(form.file_upload.data.filename)
            filepath = os.path.join(directory, filename)
            form.file_upload.data.save(filepath)

            doc = Document(str(form.file_upload.data.filename), str(form.description.data), str(filepath))
            db.session.add(doc)
            i.document.append(doc)
            db.session.commit()

            flash("Document uploaded successfully", "success")
            return redirect(url_for("investigation", name=i.investigation_name, iid=iid))

    return render_template("add_document.html", form=form, iid=iid)


@app.route("/remove_document/<int:iid>|<int:did>")
@login_required("ANY")
def remove_document(iid=-1, did=-1):
    if iid == -1 or did == -1:
        flash("Incorrect arguments for query provided", "error")
        return redirect(url_for("index"))

    from app.models import Document, Investigation
    doc = Document.query.filter_by(id=did, investigation_id=iid).first()
    path = doc.location

    if os.path.exists(path):
        os.remove(path)

    i = Investigation.query.filter_by(investigation_id=iid).first()
    db.session.delete(doc)
    db.session.commit()
    return redirect(url_for("investigation", name=i.investigation_name, iid=iid))


# TODO START
@app.route("/link_project/<int:iid>", methods=["GET", "POST"])
@login_required("ANY")
def link_project(iid=-1):
    if iid == -1:
        flash("Incorrect arguments for query provided", "error")
        return redirect(url_for("index"))

    return redirect(url_for("index"))


@app.route("/projects", methods=["GET", "POST"])
@app.route("/projects/<int:page>", methods=["GET", "POST"])
@login_required("ANY")
def projects(page=1):
    from app.models import Project
    p = Project.query.paginate(page=page, per_page=20)
    return render_template("projects.html", page=page, projects=p)


# TODO START
@app.route("/runs", methods=["GET", "POST"])
@app.route("/runs/<int:page>", methods=["GET", "POST"])
@login_required("ANY")
def runs(page=1):
    from app.models import Run
    r = Run.query.paginate(page=page, per_page=20)
    return render_template("runs.html", page=page, runs=r)


# TODO START
@app.route("/run/<int:rid>|<string:type>")
def run(rid=-1, type="none"):
    from app.models import Run

    if type == "none":
        r = Run.query.filter_by(id=rid).first()
        type = r.type

    if type == "sequencing":
        flash("Run view is not yet implemented", "warning")
        # TODO: Display sequencing run information

    elif type == "flow_cytometry":
        flash("Run view is not yet implemented", "warning")
        # TODO: Display flow cytometry run information

    flash("There was an error whilst navigating", "error")
    return render_template("index.html")


# TODO FINISH
@app.route("/input_sequencing_run", methods=["GET", "POST"])
@app.route("/input_sequencing_run/<string:type>", methods=["GET", "POST"])
@login_required("ANY")
def input_sequencing_run(type="none"):
    form = forms.UploadSequencingProjectForm(prefix="form")
    form2 = forms.NewSequencingProjectForm(prefix="form2")
    if request.method == "GET":
        return render_template("new_sequencing_project.html", form=form, form2=form2)

    elif type == "download":
        # TODO: Download the file
        data = [
            ["Sequencing Run Name", ""],
            ["Flow Cell ID", ""],
            ["Start Date (YYYY-MM-DD format)", ""],
            ["Completion Date (YYYY-MM-DD format)", ""],
            ["Genomics Lead", ""],
            ["Data Location", ""],
            ["Number of Index Tag Cycles", ""],
            ["Number of Read Cycles", ""],
            ["Paired End (y/n)", ""]
        ]
        return excel.make_response_from_array(data, "xls")

    elif type == "upload":
        if form.validate_on_submit():
            if form.file_upload.has_file():
                submission = jsonify(request.files['file_upload'])
                submission

                # TODO: Load file into memory
                # TODO: Hand off file to db model for parsing
                # TODO: Save and commit

            flash("Document uploaded successfully", "success")
            # TODO: Switch to run
            return redirect(url_for("runs", page=1))

    elif type == "manual":
        if form2.validate_on_submit():
            from app.models import SequencingRun
            s = SequencingRun()
            s.name = form2.sequence_run_name.data
            s.start_date = form2.start_date.data
            s.completion_date = form2.completion_date.data
            s.data_location = form2.data_location.data
            s.flow_cell_id = form2.flow_cell_id.data
            s.genomics_lead = form2.genomics_lead.data
            s.index_tag_cycles = form2.index_tag_cycles.data
            s.read_cycles = form2.read_cycles.data
            s.paired_end = form2.paired_end.data

            db.session.add(s)
            db.session.commit()

            flash("Sequencing run submitted successfully", "success")
            return redirect(url_for("run", rid=s.id, type=s.type))

    flash("There was an error whilst navigating.", "error")
    return render_template("new_sequencing_project.html", form=form, form2=form2)


# TODO Start
@app.route("/input_flow_cytometry_run")
@login_required("ANY")
def input_flow_cytometry_run():
    flash("New sequencing run page is still in development", "warning")
    return redirect(url_for("empty"))

# @login_required to secure
