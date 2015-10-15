import os
import flask.ext.excel as excel
import pyexcel
import pyexcel.ext.xls
import pyexcel.ext.xlsx

from app import app, db, forms, login_required
from flask import render_template, request, flash, redirect, url_for, g, jsonify
from flask.ext.login import login_user, logout_user

__author__ = 'jon'


@app.route("/")
@app.route("/index")
def index():
    if g.user is not None and g.user.is_authenticated:
        return redirect(url_for("activity"))

    return render_template("index.html", title="Home", user=g.user)


@app.route("/about")
def about():
    return render_template("index.html", title="About", user=g.user)


@app.route("/data_processing")
def data_processing():
    return render_template("index.html", title="Data Processing", user=g.user)


@app.route("/data_management")
def data_management():
    return render_template("index.html", title="Data Management", user=g.user)


@app.route("/data_monitoring")
def data_monitoring():
    return render_template("index.html", title="Data Monitoring", user=g.user)


@app.route("/data_analysis")
def data_analysis():
    return render_template("index.html", title="Data Assessment", user=g.user)


@app.route("/empty")
def empty():
    return render_template("empty.html", title="Down the rabbit hole!")


@app.route("/register", methods=["GET", "POST"])
def register():
    if g.user is not None and g.user.is_authenticated:
        return redirect(url_for("index"))

    form = forms.RegisterForm()
    if request.method == "GET":
        return render_template("register.html", title="Sign Up!", form=form)

    else:
        if not form.validate():
            return render_template("register.html", title="Sign Up!", form=form)

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
        return render_template("login.html", title="Login", form=form)

    else:
        if not form.validate_on_submit():
            return render_template("login.html", title="Login", form=form)

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


@app.route("/terms_and_conditions")
def terms_and_conditions():
    return render_template("terms_and_conditions.html", title="Terms and Conditions")


@app.route("/activity")
@login_required("ANY")
def activity():
    return render_template("activity.html", title="Interaction Panel")


@app.route("/investigations", methods=["GET", "POST"])
@app.route("/investigations/<int:page>", methods=["GET", "POST"])
@login_required("ANY")
def investigations(page=1):
    from app.models import Investigation
    i = Investigation.query.paginate(page=page, per_page=20)
    return render_template("investigations.html", title="My Investigations", page=page, obs=i)


@app.route("/new_investigation", methods=["GET", "POST"])
@login_required("ANY")
def new_investigation():
    form = forms.NewInvestigationForm()
    if request.method == "GET":
        return render_template("new_investigation.html", title="New Investigation", form=form)

    else:
        if form.validate_on_submit():
            from app.models import Investigation
            p = Investigation(str(form.investigation_name.data), str(form.investigation_lead.data))
            db.session.add(p)
            db.session.commit()
            flash("Investigation successfully registered!", "info")
            return redirect(url_for("investigations"))

        return render_template("new_investigation.html", title="New Investigation", form=form)


@app.route("/investigation/<string:name>|<int:iid>", methods=["GET", "POST"])
@login_required("ANY")
def investigation(name="", iid=-1):
    if name == "" or iid == -1:
        flash("Incorrect arguments for query provided.", "error")
        return redirect(url_for("index"))

    from app.models import Investigation
    i = Investigation.query.filter_by(investigation_id=iid, investigation_name=name).first()
    return render_template("investigation.html", title="My Investigations", investigation=i)


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

    return render_template("add_document.html", title="Add Document", form=form, iid=iid)


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


@app.route("/link_project/<int:iid>|<int:page>|<int:pid>", methods=["GET", "POST"])
@app.route("/link_project/<int:iid>|<int:page>", methods=["GET", "POST"])
@app.route("/link_project/<int:iid>", methods=["GET", "POST"])
@login_required("ANY")
def link_project(iid=-1, page=1, pid=-1):
    if pid is not -1 and iid is not -1:
        # Link the project and return
        from app.models import Investigation, Project
        i = Investigation.query.filter_by(investigation_id=iid).first()
        p = Project.query.filter_by(id=pid).first()
        i.project.append(p)
        return render_template("investigation.html", title="My Investigations", investigation=i)

    from app.models import Project
    p = Project.query.paginate(page=page, per_page=20)
    return render_template("link_project.html", title="Link Project", page=page, projects=p, iid=iid)


@app.route("/runs", methods=["GET", "POST"])
@app.route("/runs/<int:page>", methods=["GET", "POST"])
@login_required("ANY")
def runs(page=1):
    from app.models import Run
    r = Run.query.paginate(page=page, per_page=20)
    return render_template("runs.html", title="My Runs", page=page, obs=r)


# TODO FINISH
@app.route("/run/<int:rid>|<string:type>")
@login_required("ANY")
def run(rid=-1, type="none"):
    from app.models import Run

    if type == "none":
        r = Run.query.filter_by(id=rid).first()
        type = r.type

    if type == "sequencing":
        return redirect(url_for("sequencing_run", rid=rid))

    elif type == "flow_cytometry":
        flash("Run view is not yet implemented", "warning")
        # TODO: Display flow cytometry run information

    flash("There was an error whilst navigating", "error")
    return render_template("empty.html", title="Down the rabbit hole!")


@app.route("/sequencing_run/<int:rid>")
@login_required("ANY")
def sequencing_run(rid=-1):
    if rid == -1:
        flash("Incorrect arguments for query provided.", "error")
        return redirect(url_for("index"))

    from app.models import SequencingRun
    s = SequencingRun.query.filter_by(id=rid).first()
    if s is None:
        flash("Incorrect arguments for query provided.", "error")
        return redirect(url_for("index"))

    return render_template("sequencing_run.html", title="Sequencing Run", run=s)


@app.route("/input_sequencing_run", methods=["GET", "POST"])
@app.route("/input_sequencing_run/<string:type>", methods=["GET", "POST"])
@login_required("ANY")
def input_sequencing_run(type="none"):
    form = forms.UploadSequencingRunForm(prefix="form")
    form2 = forms.NewSequencingProjectForm(prefix="form2")
    if request.method == "GET":
        if type == "download":
            data = [
                ["Sequencing Run Name", ""],
                ["Flow Cell ID", ""],
                ["Start Date (YYYY-MM-DD format)", ""],
                ["Completion Date (YYYY-MM-DD format)", ""],
                ["Genomics Lead", ""],
                ["Data Location", ""],
                ["Number of Cycles for Index Tag 1", ""],
                ["Number of Cycles for Index Tag 2", ""],
                ["Number of Cycles for Read 1", ""],
                ["Number of Cycles for Read 2", ""],
                ["Paired End (Yes/No)", ""]
            ]
            response = excel.make_response_from_array(data, "xls")
            response.headers["Content-Disposition"] = "attachment; filename=template.xls"
            return response
        return render_template("input_sequencing_run.html", title="Input Sequencing Run", form=form, form2=form2)

    elif type == "upload":
        if form.validate_on_submit():
            if form.file_upload.has_file():
                filename = request.files["form-file_upload"].filename
                extension = filename.split(".")[1]

                sheet = pyexcel.load_from_memory(extension, request.files["form-file_upload"].read())
                data = dict(pyexcel.to_array(sheet))

                empty = [
                    ["Sequencing Run Name", ""],
                    ["Flow Cell ID", ""],
                    ["Start Date (YYYY-MM-DD format)", ""],
                    ["Completion Date (YYYY-MM-DD format)", ""],
                    ["Genomics Lead", ""],
                    ["Data Location", ""],
                    ["Number of Cycles for Index Tag 1", ""],
                    ["Number of Cycles for Index Tag 2", ""],
                    ["Number of Cycles for Read 1", ""],
                    ["Number of Cycles for Read 2", ""],
                    ["Paired End (Yes/No)", ""]
                ]

                for i in range(0, len(data)):
                    if empty[i][0] in data:
                        if data.get(empty[i][0]) == "":
                            flash("File was missing data for: " + empty[i][0], "error")
                            return render_template("input_sequencing_run.html", title="Input Sequencing Run", form=form,
                                                   form2=form2)

                from app.models import SequencingRun
                s = SequencingRun()
                s.name = data.get("Sequencing Run Name")
                s.start_date = data.get("Start Date (YYYY-MM-DD format)")
                s.completion_date = data.get("Completion Date (YYYY-MM-DD format)")
                s.data_location = data.get("Data Location")
                s.flow_cell_id = data.get("Flow Cell ID")
                s.genomics_lead = data.get("Genomics Lead")
                s.index_tag_cycles = data.get("Number of Cycles for Index Tag 1")
                s.index_tag_cycles_2 = data.get("Number of Cycles for Index Tag 2")
                s.read_cycles = data.get("Number of Cycles for Read 1")
                s.read_cycles_2 = data.get("Number of Cycles for Read 2")
                s.paired_end = data.get("Paired End (Yes/No)")

                db.session.add(s)
                db.session.commit()

                flash("Document uploaded successfully", "success")
                return redirect(url_for("run", rid=s.id, type=s.type))

        flash("Missing file information", "error")

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
    return render_template("input_sequencing_run.html", title="Input Sequencing Run", form=form, form2=form2)


@app.route("/projects", methods=["GET", "POST"])
@app.route("/projects/<int:page>", methods=["GET", "POST"])
@login_required("ANY")
def projects(page=1):
    from app.models import Project
    p = Project.query.paginate(page=page, per_page=20)
    return render_template("projects.html", title="My Projects", page=page, obs=p)


# TODO FINISH
@app.route("/project/<int:pid>|<string:type>")
@login_required("ANY")
def project(pid=-1, type="none"):
    from app.models import Run

    if type == "none":
        r = Run.query.filter_by(id=pid).first()
        type = r.type

    if type == "sequencing":
        return redirect(url_for("sequencing_project", pid=pid))

    elif type == "flow_cytometry":
        flash("Run view is not yet implemented", "warning")
        # TODO: Display flow cytometry run information

    flash("There was an error whilst navigating", "error")
    return render_template("empty.html", title="Down the rabbit hole!")


@app.route("/sequencing_project/<int:pid>")
@login_required("ANY")
def sequencing_project(pid=-1):
    from app.models import SequencingProject
    p = SequencingProject.query.filter_by(id=pid).first()
    return render_template("sequencing_project.html", title="Sequencing Project", project=p)


@app.route("/input_sequencing_project/<int:rid>|<string:type>", methods=["GET", "POST"])
@app.route("/input_sequencing_project", methods=["GET", "POST"])
@login_required("ANY")
def input_sequencing_project(type="", rid=-1):
    form = forms.UploadSequencingProjectForm()
    if request.method == "GET":
        if type == "download":
            # We could hold this on disk but it does not seem worth the io
            data = [
                [
                    "Internal Sample Name", "Customer Sample Name", "Sequencing Project Group Name", "Sequencing Type",
                    "Customer Name", "Lane Number", "Sequencing Concentration", "PhiXSpiked", "Spike", "Spike Ratio",
                    "Index 1 Tag Sequence", "Index 2 Tag Sequence", "Index 1 Tag ID", "Index 2 Tag ID",
                    "Index 1 Tag Kit ID", "Index 2 Tag Kit ID", "Adaptor Sequence"
                ]
            ]

            # dump the above array into an excel file for the user to download
            response = excel.make_response_from_array(data, "xls")
            response.headers["Content-Disposition"] = "attachment; filename=template.xls"
            return response
        return render_template("input_sequencing_project.html", title="Input Sequencing Project", form=form, rid=rid)

    elif type == "upload":
        if form.validate_on_submit():
            if form.file_upload.has_file():
                filename = request.files["file_upload"].filename
                extension = filename.split(".")[1]

                # Magic to get an csv/excel, transpose it and assign the first value as the dict key
                sheet = pyexcel.load_from_memory(extension, request.files["file_upload"].read())
                raw = pyexcel.to_array(sheet)
                transposed = list(zip(*raw))
                data = dict([(k[0], k[1:]) for k in transposed])

                # List of expected content - TODO: update with what information is absolutely necessary
                empty = [
                    ["Internal Sample Name", "y"],
                    ["Customer Sample Name", "y"],
                    ["Sequencing Project Group Name", "y"],
                    ["Sequencing Type", "y"],
                    ["Customer Name", "y"],
                    ["Lane Number", "y"],
                    ["Sequencing Concentration", "y"],
                    ["PhiXSpiked", "y"],
                    ["Spike", "y"],
                    ["Spike Ratio", "y"],
                    ["Index 1 Tag Sequence", "y"],
                    ["Index 2 Tag Sequence", "n"],
                    ["Index 1 Tag ID", "y"],
                    ["Index 2 Tag ID", "n"],
                    ["Index 1 Tag Kit ID", "y"],
                    ["Index 2 Tag Kit ID", "n"],
                    ["Adaptor Sequence", "y"]
                ]

                # Loop through the content, if an item is expected then throw an error when it is absent
                for i in range(0, len(data)):
                    if empty[i][1] == "y":
                        if data.get(empty[i][0]) is None:
                            flash("Could not identify the column: " + empty[i][0] + " in your submission", "error")
                            return render_template("input_sequencing_project.html", title="Input Sequencing Project",
                                                   form=form, rid=rid)

                        for j in range(0, len(data.get("Internal Sample Name"))):
                            if data.get(empty[i][0])[j] is None or data.get(empty[i][0])[j] == "":
                                flash("File was missing data for: " + empty[i][0] + " on line: " + str(
                                    j) + " this information is considered essential.", "error")
                                return render_template("input_sequencing_project.html",
                                                       title="Input Sequencing Project", form=form, rid=rid)

                # Get or create the run
                from app.models import SequencingRun
                r = SequencingRun.query.filter_by(id=rid).first()
                if r is None:
                    flash("Could not identify the parent run", "error")
                    return render_template("input_sequencing_project", title="Input Sequencing Project", form=form)

                for i in range(0, len(data.get("Internal Sample Name"))):

                    # Create a sample which must be unique in the context of the run, this is done first so we can error
                    # out and leave the db intact if something is wrong
                    from app.models import SequencingSample, Tag
                    s = r.sample.filter_by(internal_sample_name=data.get("Internal Sample Name")[i]).first()
                    if s is None:
                        s = SequencingSample()
                        s.internal_sample_name = data.get("Internal Sample Name")[i]
                        s.customer_sample_name = data.get("Customer Sample Name")[i]
                        s.adaptor_sequence = data.get("Adaptor Sequence")[i]

                        # Not all samples use both tags
                        s.index_tag.append(Tag(data.get("Index 1 Tag ID")[i], data.get("Index 1 Tag Kit ID")[i],
                                               data.get("Index 1 Tag Sequence")[i], True))
                        if data.get("Index 2 Tag ID")[i] is not None and data.get("Index 2 Tag ID") is not "":
                            s.index_tag.append(Tag(data.get("Index 2 Tag ID")[i], data.get("Index 2 Tag Kit ID")[i],
                                                   data.get("Index 2 Tag Sequence")[i], False))

                        db.session.add(s)
                        db.session.commit()

                    # if it is already present we error out of the submission
                    else:
                        flash(
                            "Your submission contained a sample name that already exists within the context of this" +
                            " run. The internal sample name in question was: " +
                            str(data.get("Internal Sample Name")[i]) +
                            " please re-submit when corrections are made", "error")
                        return render_template("input_sequencing_project.html", title="Input Sequencing Project",
                                               form=form, rid=rid)

                    r.sample.append(s)

                    # Create a customer if their name has not been used before, otherwise we will blindly update the
                    # existing customer's information - not a lot we can do here
                    from app.models import Customer
                    c = Customer.query.filter_by(name=data.get("Customer Name")[i]).first()
                    if c is None:
                        c = Customer(data.get("Customer Name")[i])

                        db.session.add(c)
                        db.session.commit()

                    c.sample.append(s)  # We can be certain of this because if it is not new we error out

                    # Create a project grouping, the name of which must be unique within the context of the run
                    # No handling for repeats as we update the existing record if the project group is already present
                    from app.models import SequencingProject
                    p = r.project.filter_by(name=data.get("Sequencing Project Group Name")[i]).first()
                    if p is None:
                        p = SequencingProject()
                        p.name = data.get("Sequencing Project Group Name")[i]
                        p.sequencing_type = data.get("Sequencing Type")[i]

                        db.session.add(p)
                        db.session.commit()

                    r.project.append(p)
                    c.project.append(p)
                    p.sample.append(s)

                    # Create or retrieve lane information which must be unique in the context of the run
                    from app.models import Lane
                    l = r.lane.filter_by(number=data.get("Lane Number")[i]).first()
                    if l is None:
                        l = Lane(data.get("Lane Number")[i])
                        l.set_sequencing_concentration(data.get("Sequencing Concentration")[i])
                        l.phi_x_spiked = data.get("PhiXSpiked")[i]
                        l.spike = data.get("Spike")[i]
                        l.spikeRation = data.get("Spike Ratio")[i]

                        db.session.add(l)
                        db.session.commit()

                    r.lane.append(l)
                    s.lane = l  # This is safe because we know the sample is new

                # Save the run (should append any nested creations etc) - this probably isn't necessary but I cannot be
                # sure right now!
                db.session.add(r)
                db.session.commit()

                flash("Document uploaded successfully", "success")
                return redirect(url_for("sequencing_run", rid=r.id, type=""))

        # Handle invalid form submissions
        flash("Missing file information", "error")

    flash("There was an error whilst navigating.", "error")
    return render_template("input_sequencing_project.html", title="Input Sequencing Project", form=form, rid=rid)


# TODO Start
@app.route("/input_flow_cytometry_run")
@login_required("ANY")
def input_flow_cytometry_run():
    flash("New flow cytometry run page is still in development", "warning")
    return redirect(url_for("empty"))


@app.route("/demultiplex/<int:rid>|<int:pid>", methods=["GET", "POST"])
@app.route("/demultiplex", methods=["GET", "POST"])
@login_required("ANY")
def demultiplex(rid=-1, pid=-1):
    flash("demultiplex information is not yet implemented", "warning")
    return redirect(url_for("sequencing_run", rid=rid))


@app.route("/fast_qc/<int:rid>|<int:pid>", methods=["GET", "POST"])
@app.route("/fast_qc", methods=["GET", "POST"])
@login_required("ANY")
def fast_qc(rid=-1, pid=-1):
    flash("fast qc information is not yet implemented", "warning")
    return redirect(url_for("sequencing_run", rid=rid))


@app.route("/post_align/<int:rid>|<int:pid>", methods=["GET", "POST"])
@app.route("/post_align", methods=["GET", "POST"])
@login_required("ANY")
def post_align(rid=-1, pid=-1):
    flash("post align is not yet implemented", "warning")
    return redirect(url_for("sequencing_run", rid=rid))


@app.route("/sftp/<int:rid>|<int:pid>", methods=["GET", "POST"])
@app.route("/sftp", methods=["GET", "POST"])
@login_required("ANY")
def sftp(rid=-1, pid=-1):
    flash("sftp not yet implemented", "warning")
    return redirect(url_for("sequencing_run", rid=rid))

# @login_required to secure
