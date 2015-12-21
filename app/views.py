import os
import subprocess

import config
import flask.ext.excel as excel
import pyexcel
import pyexcel.ext.xls
import pyexcel.ext.xlsx

from app import app, db, forms, login_required, utils, login_manager
from flask import render_template, request, flash, redirect, url_for, g
from flask.ext.login import login_user, logout_user, current_user

__author__ = 'jon'


@app.route("/")
@app.route("/index")
def index():
    if g.user is not None and g.user.is_authenticated:
        return redirect(url_for("activity"))

    return render_template("welcome.html", title="Home", user=g.user)


@app.route("/about")
def about():
    return render_template("welcome.html", title="About", user=g.user)


@app.route("/data_processing")
def data_processing():
    return render_template("data_processing.html", title="Data Processing", user=g.user)


@app.route("/data_management")
def data_management():
    return render_template("data_management.html", title="Data Management", user=g.user)


@app.route("/data_monitoring")
def data_monitoring():
    return render_template("data_monitoring.html", title="Data Monitoring", user=g.user)


@app.route("/data_analysis")
def data_analysis():
    return render_template("data_analysis.html", title="Data Assessment", user=g.user)


@app.route("/empty")
def empty():
    return render_template("empty.html", title="Down the rabbit hole!")


@app.route("/login", methods=["GET", "POST"])
def login():
    if g.user is not None and g.user.is_authenticated:
        return redirect(url_for("index"))

    form = forms.LoginForm()
    if request.method == "GET":
        return render_template("login.html", title="Login", form=form)

    else:
        if not form.validate_on_submit():
            utils.flash_errors(form)
            return render_template("login.html", title="Login", form=form)

        else:
            from app.models import User
            user = User.query.filter_by(login_name=str(form.username.data)).first()
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


@app.route("/administrate")
@login_required("Site Admin", "Group Admin")
def administrate():
    return render_template("administrate_welcome.html", title="Administrate Panel")


@app.route("/refresh_pipelines")
@login_required("Site Admin")
def refresh_pipelines():
    utils.refresh_pipelines()
    return redirect(url_for("administrate"))


@app.route("/show_groups")
@app.route("/show_groups/<int:page>")
@login_required("Site Admin")
def show_groups(page=1):
    from app.models import Group
    g = Group.query.paginate(page=page, per_page=20)
    return render_template("groups.html", title="Groups", page=page, obs=g)


@app.route("/show_users")
@app.route("/show_users/<int:page>")
@login_required("Site Admin", "Group Admin")
def show_users(page=1):
    from app.models import User
    if current_user.get_role() == "Site Admin":
        u = User.query.paginate(page=page, per_page=20)

    else:
        g = current_user.group
        u = User.query.filter_by(group_id=g.id).paginate(page=page, per_page=20)

    return render_template("users.html", title="Users", page=page, obs=u)


@app.route("/show_customers")
@app.route("/show_customers/<int:page>")
@login_required("Site Admin", "Group Admin")
def show_customers(page=1):
    from app.models import Customer
    if current_user.get_role() == "Site Admin":
        c = Customer.query.paginate(page=page, per_page=20)

    else:
        g = current_user.group
        c = Customer.query.filter_by(group_id=g.id).paginate(page=page, per_page=20)

    return render_template("users.html", title="Users", page=page, obs=c)


@app.route("/add_group", methods=["GET", "POST"])
@login_required("Site Admin")
def add_group():
    form = forms.CreateGroupForm()
    if request.method == "GET":
        return render_template("add_group.html", title="Group Creation", form=form)

    else:
        if not form.validate():
            utils.flash_errors(form)
            return render_template("add_group.html", title="Group Creation", form=form)

        else:
            from app.models import Group, User

            admin = User()
            admin.login_name = str(form.admin_login.data)
            admin.set_password(str(form.admin_password.data))
            admin.email = str(form.admin_email.data)
            admin.set_role("Group Admin")

            group = Group()
            group.name = str(form.group_name.data)
            group.member.append(admin)

            db.session.add(admin)
            db.session.add(group)
            db.session.commit()
            return redirect(url_for("show_groups"))


@app.route("/add_user", methods=["GET", "POST"])
@login_required("Group Admin")
def add_user():
    form = forms.CreatePerson()
    if request.method == "GET":
        return render_template("add_person.html", title="Add User", type="User", form=form)

    else:
        if not form.validate():
            utils.flash_errors(form)
            return render_template("add_person.html", title="Add User", type="User", form=form)

        else:
            from app.models import User
            user = User()
            user.login_name = str(form.login_name.data)
            user.set_password(str(form.login_password.data))
            user.email = str(form.login_email.data)

            current_user.group.member.append(user)
            db.session.add(user)
            db.session.add(current_user.group)
            db.session.commit()
            return redirect(url_for("show_users"))


@app.route("/add_customer", methods=["GET", "POST"])
@login_required("Group Admin")
def add_customer():
    form = forms.CreatePerson()
    if request.method == "GET":
        return render_template("add_person.html", title="Add User", type="Customer", form=form)

    else:
        if not form.validate():
            utils.flash_errors(form)
            return render_template("add_person.html", title="Add User", type="Customer", form=form)

        else:
            from app.models import Customer
            user = Customer()
            user.login_name = str(form.login_name.data)
            user.set_password(str(form.login_password.data))
            user.email = str(form.login_email.data)

            current_user.group.member.append(user)
            db.session.add(user)
            db.session.add(current_user.group)
            db.session.commit()
            return redirect(url_for("show_customers"))


@app.route("/investigations", methods=["GET", "POST"])
@app.route("/investigations/<int:page>", methods=["GET", "POST"])
@login_required("ANY")
def investigations(page=1):
    from app.models import Investigation
    if current_user.get_role() == "Site Admin":
        i = Investigation.query.paginate(page=page, per_page=20)

    else:
        i = utils.get_allowed_investigations_query().paginate(page=page, per_page=20)

    return render_template("investigations.html", title="My Investigations", page=page, obs=i)


@app.route("/new_investigation", methods=["GET", "POST"])
@login_required("ANY")
def new_investigation():
    form = forms.NewInvestigationForm()
    if request.method == "GET":
        return render_template("new_investigation.html", title="New Investigation", form=form)

    else:
        if form.validate_on_submit():
            utils.create_investigation(str(form.investigation_name.data), str(form.investigation_lead.data))
            flash("Investigation successfully registered!", "info")
            return redirect(url_for("investigations"))

        return render_template("new_investigation.html", title="New Investigation", form=form)


@app.route("/investigation/<iid>", methods=["GET", "POST"])
@login_required("ANY")
def investigation(iid=""):
    if iid == "":
        flash("Incorrect arguments for query provided.", "error")
        return redirect(url_for("index"))

    i = utils.get_allowed_investigation_by_display_key(iid)
    return render_template("investigation.html", title="My Investigations", investigation=i)


@app.route("/add_document/<iid>", methods=["GET", "POST"])
@login_required("ANY")
def add_document(iid=""):
    if iid == "":
        flash("Incorrect arguments for query provided", "error")
        return redirect(url_for("index"))

    form = forms.AddDocumentForm()
    if request.method == "POST":
        if form.validate_on_submit():
            # Handle no investigation to link to
            i = utils.get_allowed_investigation_by_display_key(iid)
            if not i:
                flash("Attempted to add a document to an investigation without an investigation parent.", "error")
                return redirect(url_for("investigations"))

            # Handle a bad directory structure - should never happen
            directory = i.validate_investigation_directory()
            if directory is None:
                flash("Could not create investigation directory.", "error")
                return redirect(url_for("investigations"))

            # Handle maliciously named files (i.e. ../..)
            from werkzeug.utils import secure_filename
            filename = secure_filename(form.file_upload.data.filename)
            filepath = os.path.join(directory, filename)

            # Handle a document already existing
            if os.path.exists(filepath):
                flash("A document with this location (i.e. filename) already exists", "error")
                return redirect(url_for("investigations"))

            # Save the file to the given path
            form.file_upload.data.save(filepath)

            # Save the document to the db
            utils.create_document(iid, str(form.file_upload.data.filename), str(form.description.data), str(filepath))

            # Inform and redirect
            flash("Document uploaded successfully", "success")
            return redirect(url_for("investigation", iid=iid))

    # Fail scenario
    return render_template("add_document.html", title="Add Document", form=form, iid=iid)


@app.route("/remove_document/<iid>|<did>")
@login_required("ANY")
def remove_document(iid="", did=""):
    if iid == "" or did == "":
        flash("Incorrect arguments for query provided", "error")
        return redirect(url_for("index"))

    doc = utils.get_allowed_document_by_display_key(did)
    if doc.investigation_id is not did:
        flash("Incorrect arguments for query provided", "error")
        return redirect(url_for("index"))

    path = doc.location
    if os.path.exists(path):
        os.remove(path)

    utils.remove_document(did)
    return redirect(url_for("investigation", iid=iid))


@app.route("/link_sample_group/<iid>|<gid>|<int:page>", methods=["GET", "POST"])
@app.route("/link_sample_group/<iid>|<int:page>", methods=["GET", "POST"])
@app.route("/link_sample_group/<iid>", methods=["GET", "POST"])
@login_required("ANY")
def link_sample_group(iid="", gid="", page=1):
    if gid is not "" and iid is not "":
        # Link the project and return
        i = utils.get_allowed_investigation_by_display_key(iid)
        if utils.link_sample_group_to_investigation(iid, gid):
            return render_template("investigation.html", title="My Investigations", investigation=i)

        else:
            flash("Cannot find or you do not have access to the provided investigation and/or sample group", "error")
            return render_template("investigation.html", title="My Investigations", investigation=i)

    p = utils.get_allowed_sample_groups_query().paginate(page=page, per_page=20)
    return render_template("link_sample_group.html", title="Link Sample Group", page=page, obs=p, iid=iid)


@app.route("/submissions", methods=["GET", "POST"])
@app.route("/submissions/<int:page>", methods=["GET", "POST"])
@login_required("ANY")
def submissions(page=1):
    s = utils.get_allowed_submissions_query().paginate(page=page, per_page=20)
    return render_template("submissions.html", title="My Runs", page=page, obs=s)


# TODO FINISH
@app.route("/submission/<sid>|<string:type>")
@login_required("ANY")
def submission(sid="", type="none"):
    if type == "none":
        s = utils.get_allowed_submission_by_display_key(sid)
        type = s.type

    if type == "Sequencing":
        return redirect(url_for("sequencing_submission", sid=sid))

    elif type == "Flow Cytometry":
        flash("Submission view is not yet implemented", "warning")
        # TODO: Display flow cytometry run information

    flash("There was an error whilst navigating", "error")
    return render_template("empty.html", title="Down the rabbit hole!")


@app.route("/sequencing_submission/<sid>")
@login_required("ANY")
def sequencing_submission(sid=""):
    if sid == "":
        flash("Incorrect arguments for query provided.", "error")
        return redirect(url_for("index"))

    s = utils.get_allowed_submission_by_display_key(sid)
    if s is None:
        flash("Incorrect arguments for query provided.", "error")
        return redirect(url_for("index"))

    return render_template("sequencing_submission.html", title="Sequencing Submission", submission=s)


@app.route("/input_sequencing_submission", methods=["GET", "POST"])
@app.route("/input_sequencing_submission/<string:type>", methods=["GET", "POST"])
@login_required("ANY")
def input_sequencing_submission(type="none"):
    if current_user.type == "Customer":
        flash("Customers may not access this page", "error")
        return login_manager.unauthorized()

    form = forms.UploadSequencingSubmissionForm(prefix="form")
    form2 = forms.NewSequencingSubmissionForm(prefix="form2")
    from app.static.io import sequencing_submission_template
    if request.method == "GET":
        if type == "download":
            response = excel.make_response_from_book_dict(sequencing_submission_template.submission_data, "xls")
            response.headers["Content-Disposition"] = "attachment; filename=template_for_sequencing_submission.xls"
            return response
        return render_template("input_sequencing_submission.html", title="Input Sequencing Submission", form=form, form2=form2)

    elif type == "upload":
        if form.validate_on_submit():
            if form.file_upload.has_file():
                filename = request.files["form-file_upload"].filename
                extension = filename.split(".")[1]

                try:
                    sheet = pyexcel.load_from_memory(extension, request.files["form-file_upload"].read(), sheetname="Submission Data")

                except ValueError:
                    flash("Could not find the Submission Data sheet in the provided excel file", "error")
                    return render_template("input_sequencing_submission.html", title="Input Sequencing Submission", form=form, form2=form2)

                array = pyexcel.to_array(sheet)
                for entry in array:
                    del entry[1]

                if len(array[0]) is not 2:
                    flash("Either no information was provided, or too much. We are only expecting one column of data beyond our conventional information (i.e. 3rd column is the last column!", "error")
                    return render_template("input_sequencing_submission.html", title="Input Sequencing Submission", form=form, form2=form2)

                data = dict(array)

                error = sequencing_submission_template.validate_data_sheet(data)
                if error is not None:
                    flash("File was missing data for: " + error, "error")
                    return render_template("input_sequencing_submission.html", title="Input Sequencing Submission", form=form,
                                           form2=form2)

                s = sequencing_submission_template.build_submission_entry(data)
                if s is None:
                    return render_template("input_sequencing_submission.html", title="Input Sequencing Submission", form=form,
                                           form2=form2)

                flash("Document uploaded successfully", "success")
                return redirect(url_for("submission", sid=s.display_key, type=s.type))

            else:
                flash("Missing file information", "error")

        else:
            flash("Missing file information", "error")

    elif type == "manual":
        if form2.validate_on_submit():
            s = utils.create_sequencing_submission(
                form2.sequence_submission_name.data,
                form2.start_date.data,
                form2.completion_date.data,
                form2.data_location.data,
                form2.flow_cell_id.data,
                form2.genomics_lead.data,
                form2.index_tag_cycles.data,
                form2.read_cycles.data,
                form2.paired_end.data
            )

            flash("Sequencing submission submitted successfully", "success")
            return redirect(url_for("submission", sid=s.display_key, type=s.type))

    return render_template("input_sequencing_submission.html", title="Input Sequencing Submission", form=form, form2=form2)


@app.route("/sample_groups", methods=["GET", "POST"])
@app.route("/sample_groups/<int:page>", methods=["GET", "POST"])
@login_required("ANY")
def sample_groups(page=1):
    g = utils.get_allowed_sample_groups_query().paginate(page=page, per_page=20)
    return render_template("sample_groups.html", title="My Sample Groups", page=page, obs=g)


# TODO FINISH
@app.route("/sample_group/<gid>|<string:type>")
@login_required("ANY")
def sample_group(gid="", type="none"):
    if type == "none":
        p = utils.get_allowed_sample_group_by_display_key(gid)
        type = p.type

    if type == "Sequencing":
        return redirect(url_for("sequencing_sample_group", gid=gid))

    elif type == "Flow Cytometry":
        flash("Submission view is not yet implemented", "warning")
        # TODO: Display flow cytometry run information

    flash("There was an error whilst navigating", "error")
    return render_template("empty.html", title="Down the rabbit hole!")


@app.route("/sequencing_sample_group/<gid>")
@login_required("ANY")
def sequencing_sample_group(gid=""):
    g = utils.get_allowed_sample_group_by_display_key(gid)
    p = utils.get_pipelines()
    if g is None:
        flash("Invalid sample group", "warning")
        return redirect(url_for("index"))

    return render_template("sequencing_sample_group.html", title="Sequencing Sample Group", group=g, pipelines=p)


@app.route("/link_sample_group_to_client/<gid>")
@login_required("ANY")
def link_sample_group_to_client(gid=""):
    flash("Sample group to client linking not implemented yet.", "warning")
    return redirect(url_for("sequencing_sample_group"), gid)


@app.route("/input_sequencing_sample_mappings/<sid>|<string:type>", methods=["GET", "POST"])
@app.route("/input_sequencing_sample_mappings", methods=["GET", "POST"])
@login_required("ANY")
def input_sequencing_sample_mappings(type="", sid=""):
    if current_user.type == "Customer":
        flash("Customers may not access this page", "error")
        return login_manager.unauthorized()

    form = forms.UploadSequencingSampleMappingsForm()
    from app.static.io import sequencing_sample_mappings_template
    if request.method == "GET":
        if type == "download":
            response = excel.make_response_from_book_dict(sequencing_sample_mappings_template.sample_mappings_data, "xls")
            response.headers["Content-Disposition"] = "attachment; filename=template_for_sequencing_sample_mapping.xls"
            return response
        return render_template("input_sequencing_sample_mappings.html", title="Input Sequencing Sample Mappings", form=form, sid=sid)

    elif type == "upload":
        if form.validate_on_submit():
            if form.file_upload.has_file():
                filename = request.files["file_upload"].filename
                extension = filename.split(".")[1]

                # Magic to get an csv/excel, transpose it and assign the first value as the dict key
                try:
                    sheet = pyexcel.load_from_memory(extension, request.files["file_upload"].read(), "Sample Mappings Data")

                except ValueError:
                    flash("Could not find the Sample Mappings Data sheet in the provided excel file", "error")
                    return render_template("input_sequencing_sample_mappings.html", title="Input Sequencing Sample Mappings", form=form, sid=sid)

                raw = pyexcel.to_array(sheet)
                transposed = list(zip(*raw))
                data = dict([(k[0], k[2:]) for k in transposed])

                passed = sequencing_sample_mappings_template.validate_data_sheet(data)
                if not passed:
                    return render_template("input_sequencing_sample_mappings.html", title="Input Sequencing Sample Mappings", form=form, sid=sid)

                passed = sequencing_sample_mappings_template.build_sample_mappings(sid, data)
                if passed:
                    flash("Document uploaded successfully", "success")
                    return redirect(url_for("sequencing_submission", sid=sid, type=""))

                else:
                    return render_template("input_sequencing_sample_mappings.html", title="Input Sequencing Sample Mappings", form=form, sid=sid)

            else:
                flash("Missing file information", "error")

        else:
            flash("Missing file information", "error")

    return render_template("input_sequencing_sample_mappings.html", title="Input Sequencing Sample Mappings", form=form, sid=sid)


@app.route("/run_pipeline/<sample_group>|<pipeline>")
@login_required("ANY")
def run_pipeline(sample_group="", pipeline=""):
    g = utils.get_allowed_sample_group_by_display_key(sample_group)
    p = utils.get_pipeline_by_display_key(pipeline)

    i = utils.create_pipeline_instance(p)
    g.current_pipeline = i

    db.session.add(g)
    db.session.commit()

    # TODO: Start the pipeline!

    return redirect(url_for("sequencing_sample_group", gid=g.display_key))


@app.route("/pipeline/<pid>")
def pipeline(pid=""):
    p = utils.get_pipeline_by_display_key(pid)
    return redirect(url_for("empty"))


@app.route("/test_pipeline")
@login_required("Site Admin")
def test_pipeline(pid=""):
    # TODO REMOVE THIS:
    from app.models import Pipeline
    p = Pipeline.query.filter_by(name="test_pipeline").first()
    pid = p.display_key

    # Get the relevant pipeline and build an instance
    p = Pipeline.query.filter_by(display_key=pid).first()
    pi = utils.create_pipeline_instance(p)

    # Get the first module and build an instance
    m = p.module.filter_by(execution_index=0).first()
    mi = utils.create_module_instance(m, pi)

    # TODO: Generate module options, based on behaviours and query the server about it

    # Setup and Acquire Ticket
    # from app.models import Ticket
    # t = utils.create_ticket(mi)

    # Setup working directory
    # TODO: How does this behaviour change based on the plugin type

    # Setup samples csv
    # TODO: How does this behaviour change based on the csv

    # Step 4) Submit job
    shell_path = os.path.join(os.path.dirname(__file__), "static")
    shell_path = os.path.join(shell_path, "pipelines")
    shell_path = os.path.join(shell_path, "submit_job.sh")
    pipeline_path = m.executor
    process = subprocess.Popen([shell_path, "-t=A_Ticket", "-j=A_JOB", "-s=" + pipeline_path, "-w=rand", "-i=csv", "-v='a=1,b=eleven'"]).stdout

    return redirect(url_for("empty"))


# TODO Start
@app.route("/input_flow_cytometry_submission")
@login_required("ANY")
def input_flow_cytometry_submission():
    flash("New flow cytometry submission page is still in development", "warning")
    return redirect(url_for("empty"))


@app.route("/message", methods=["GET", "POST"])
def message():
    global tmp_message
    if request.method == "GET":
        m = tmp_message
        return render_template("message.html", dictionary=m)

    else:
        tmp_message = request.form
        return "<html></html>"

