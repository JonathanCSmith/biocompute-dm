import os
import subprocess
import config
import flask
import flask.ext.excel as excel
import pyexcel
import pyexcel.ext.xls
import pyexcel.ext.xlsx

from biocomputedm import forms, utils
from biocomputedm.decorators import login_required
from biocomputedm.extensions import login_manager, db
from flask import render_template, request, flash, redirect, url_for, g
from flask.ext.login import login_user, logout_user, current_user

__author__ = 'jon'

#
#

#
#
# @app.route("/investigations", methods=["GET", "POST"])
# @app.route("/investigations/<int:page>", methods=["GET", "POST"])
# @login_required("ANY")
# def investigations(page=1):
#     from biocomputedm.models import Investigation
#     if current_user.get_role() == "Site Admin":
#         i = Investigation.query.paginate(page=page, per_page=20)
#
#     else:
#         i = utils.get_allowed_investigations_query().paginate(page=page, per_page=20)
#
#     return render_template("investigations.html", title="My Investigations", page=page, obs=i)
#
#
# @app.route("/new_investigation", methods=["GET", "POST"])
# @login_required("ANY")
# def new_investigation():
#     form = forms.NewInvestigationForm()
#     if request.method == "GET":
#         return render_template("new_investigation.html", title="New Investigation", form=form)
#
#     else:
#         if form.validate_on_submit():
#             utils.create_investigation(str(form.investigation_name.data), str(form.investigation_lead.data))
#             flash("Investigation successfully registered!", "info")
#             return redirect(url_for("investigations"))
#
#         return render_template("new_investigation.html", title="New Investigation", form=form)
#
#
# @app.route("/investigation/<iid>", methods=["GET", "POST"])
# @login_required("ANY")
# def investigation(iid=""):
#     if iid == "":
#         flash("Incorrect arguments for query provided.", "error")
#         return redirect(url_for("index"))
#
#     i = utils.get_allowed_investigation_by_display_key(iid)
#     return render_template("investigation.html", title="My Investigations", investigation=i)
#
#
# @app.route("/add_document/<iid>", methods=["GET", "POST"])
# @login_required("ANY")
# def add_document(iid=""):
#     if iid == "":
#         flash("Incorrect arguments for query provided", "error")
#         return redirect(url_for("index"))
#
#     form = forms.AddDocumentForm()
#     if request.method == "POST":
#         if form.validate_on_submit():
#             # Handle no investigation to link to
#             i = utils.get_allowed_investigation_by_display_key(iid)
#             if not i:
#                 flash("Attempted to add a document to an investigation without an investigation parent.", "error")
#                 return redirect(url_for("investigations"))
#
#             # Handle a bad directory structure - should never happen
#             directory = i.validate_investigation_directory()
#             if directory is None:
#                 flash("Could not create investigation directory.", "error")
#                 return redirect(url_for("investigations"))
#
#             # Handle maliciously named files (i.e. ../..)
#             from werkzeug.utils import secure_filename
#             filename = secure_filename(form.file_upload.data.filename)
#             filepath = os.path.join(directory, filename)
#
#             # Handle a document already existing
#             if os.path.exists(filepath):
#                 flash("A document with this location (i.e. filename) already exists", "error")
#                 return redirect(url_for("investigations"))
#
#             # Save the file to the given path
#             form.file_upload.data.save(filepath)
#
#             # Save the document to the db
#             utils.create_document(iid, str(form.file_upload.data.filename), str(form.description.data), str(filepath))
#
#             # Inform and redirect
#             flash("Document uploaded successfully", "success")
#             return redirect(url_for("investigation", iid=iid))
#
#     # Fail scenario
#     return render_template("add_document.html", title="Add Document", form=form, iid=iid)
#
#
# @app.route("/remove_document/<iid>|<did>")
# @login_required("ANY")
# def remove_document(iid="", did=""):
#     if iid == "" or did == "":
#         flash("Incorrect arguments for query provided", "error")
#         return redirect(url_for("index"))
#
#     doc = utils.get_allowed_document_by_display_key(did)
#     if doc.investigation_id is not did:
#         flash("Incorrect arguments for query provided", "error")
#         return redirect(url_for("index"))
#
#     path = doc.location
#     if os.path.exists(path):
#         os.remove(path)
#
#     utils.remove_document(did)
#     return redirect(url_for("investigation", iid=iid))
#
#
# @app.route("/link_sample_group/<iid>|<gid>|<int:page>", methods=["GET", "POST"])
# @app.route("/link_sample_group/<iid>|<int:page>", methods=["GET", "POST"])
# @app.route("/link_sample_group/<iid>", methods=["GET", "POST"])
# @login_required("ANY")
# def link_sample_group(iid="", gid="", page=1):
#     if gid is not "" and iid is not "":
#         # Link the project and return
#         i = utils.get_allowed_investigation_by_display_key(iid)
#         if utils.link_sample_group_to_investigation(iid, gid):
#             return render_template("investigation.html", title="My Investigations", investigation=i)
#
#         else:
#             flash("Cannot find or you do not have access to the provided investigation and/or sample group", "error")
#             return render_template("investigation.html", title="My Investigations", investigation=i)
#
#     p = utils.get_allowed_sample_groups_query().paginate(page=page, per_page=20)
#     return render_template("link_sample_group.html", title="Link Sample Group", page=page, obs=p, iid=iid)
#
#

#
#
# # TODO FINISH
# @app.route("/submission/<sid>|<string:type>")
# @login_required("ANY")
# def submission(sid="", type="none"):
#     if type == "none":
#         s = utils.get_allowed_submission_by_display_key(sid)
#         type = s.type
#
#     if type == "Sequencing":
#         return redirect(url_for("sequencing_submission", sid=sid))
#
#     elif type == "Flow Cytometry":
#         flash("Submission view is not yet implemented", "warning")
#         # TODO: Display flow cytometry run information
#
#     flash("There was an error whilst navigating", "error")
#     return render_template("empty.html", title="Down the rabbit hole!")
#
#
# @app.route("/sequencing_submission/<sid>")
# @login_required("ANY")
# def sequencing_submission(sid=""):
#     if sid == "":
#         flash("Incorrect arguments for query provided.", "error")
#         return redirect(url_for("index"))
#
#     s = utils.get_allowed_submission_by_display_key(sid)
#     if s is None:
#         flash("Incorrect arguments for query provided.", "error")
#         return redirect(url_for("index"))
#
#     return render_template("sequencing_submission.html", title="Sequencing Submission", submission=s)
#
#
# @app.route("/input_sequencing_submission", methods=["GET", "POST"])
# @app.route("/input_sequencing_submission/<string:type>", methods=["GET", "POST"])
# @login_required("ANY")
# def input_sequencing_submission(type="none"):
#     if current_user.type == "Customer":
#         flash("Customers may not access this page", "error")
#         return login_manager.unauthorized()
#
#     form = forms.UploadSequencingSubmissionForm(prefix="form")
#     form2 = forms.NewSequencingSubmissionForm(prefix="form2")
#     from biocomputedm.static.io import sequencing_submission_template
#     if request.method == "GET":
#         if type == "download":
#             response = excel.make_response_from_book_dict(sequencing_submission_template.submission_data, "xls")
#             response.headers["Content-Disposition"] = "attachment; filename=template_for_sequencing_submission.xls"
#             return response
#         return render_template("input_sequencing_submission.html", title="Input Sequencing Submission", form=form,
#                                form2=form2)
#
#     elif type == "upload":
#         if form.validate_on_submit():
#             if form.file_upload.has_file():
#                 filename = request.files["form-file_upload"].filename
#                 extension = filename.split(".")[1]
#
#                 try:
#                     sheet = pyexcel.load_from_memory(extension, request.files["form-file_upload"].read(),
#                                                      sheetname="Submission Data")
#
#                 except ValueError:
#                     flash("Could not find the Submission Data sheet in the provided excel file", "error")
#                     return render_template("input_sequencing_submission.html", title="Input Sequencing Submission",
#                                            form=form, form2=form2)
#
#                 array = pyexcel.to_array(sheet)
#                 for entry in array:
#                     del entry[1]
#
#                 if len(array[0]) is not 2:
#                     flash(
#                             "Either no information was provided, or too much. We are only expecting one column of data beyond our conventional information (i.e. 3rd column is the last column!",
#                             "error")
#                     return render_template("input_sequencing_submission.html", title="Input Sequencing Submission",
#                                            form=form, form2=form2)
#
#                 data = dict(array)
#
#                 error = sequencing_submission_template.validate_data_sheet(data)
#                 if error is not None:
#                     flash("File was missing data for: " + error, "error")
#                     return render_template("input_sequencing_submission.html", title="Input Sequencing Submission",
#                                            form=form,
#                                            form2=form2)
#
#                 s = sequencing_submission_template.build_submission_entry(data)
#                 if s is None:
#                     return render_template("input_sequencing_submission.html", title="Input Sequencing Submission",
#                                            form=form,
#                                            form2=form2)
#
#                 flash("Document uploaded successfully", "success")
#                 return redirect(url_for("submission", sid=s.display_key, type=s.type))
#
#             else:
#                 flash("Missing file information", "error")
#
#         else:
#             flash("Missing file information", "error")
#
#     elif type == "manual":
#         if form2.validate_on_submit():
#             s = utils.create_sequencing_submission(
#                     form2.sequence_submission_name.data,
#                     form2.start_date.data,
#                     form2.completion_date.data,
#                     form2.data_location.data,
#                     form2.flow_cell_id.data,
#                     form2.genomics_lead.data,
#                     form2.index_tag_cycles.data,
#                     form2.read_cycles.data,
#                     form2.paired_end.data
#             )
#
#             flash("Sequencing submission submitted successfully", "success")
#             return redirect(url_for("submission", sid=s.display_key, type=s.type))
#
#     return render_template("input_sequencing_submission.html", title="Input Sequencing Submission", form=form,
#                            form2=form2)
#
#
# @app.route("/sample_groups", methods=["GET", "POST"])
# @app.route("/sample_groups/<int:page>", methods=["GET", "POST"])
# @login_required("ANY")
# def sample_groups(page=1):
#     g = utils.get_allowed_sample_groups_query().paginate(page=page, per_page=20)
#     return render_template("sample_groups.html", title="My Sample Groups", page=page, obs=g)
#
#
# # TODO FINISH
# @app.route("/sample_group/<gid>|<string:type>")
# @login_required("ANY")
# def sample_group(gid="", type="none"):
#     if type == "none":
#         p = utils.get_allowed_sample_group_by_display_key(gid)
#         type = p.type
#
#     if type == "Sequencing":
#         return redirect(url_for("sequencing_sample_group", gid=gid))
#
#     elif type == "Flow Cytometry":
#         flash("Submission view is not yet implemented", "warning")
#         # TODO: Display flow cytometry run information
#
#     flash("There was an error whilst navigating", "error")
#     return render_template("empty.html", title="Down the rabbit hole!")
#
#
# @app.route("/sequencing_sample_group/<gid>")
# @login_required("ANY")
# def sequencing_sample_group(gid=""):
#     g = utils.get_allowed_sample_group_by_display_key(gid)
#     p = utils.get_pipelines()
#     if g is None:
#         flash("Invalid sample group", "warning")
#         return redirect(url_for("index"))
#
#     return render_template("sequencing_sample_group.html", title="Sequencing Sample Group", group=g, pipelines=p)
#
#
# @app.route("/link_sample_group_to_client/<gid>")
# @login_required("ANY")
# def link_sample_group_to_client(gid=""):
#     flash("Sample group to client linking not implemented yet.", "warning")
#     return redirect(url_for("sequencing_sample_group"), gid)
#
#
# @app.route("/input_sequencing_sample_mappings/<sid>|<string:type>", methods=["GET", "POST"])
# @app.route("/input_sequencing_sample_mappings", methods=["GET", "POST"])
# @login_required("ANY")
# def input_sequencing_sample_mappings(type="", sid=""):
#     if current_user.type == "Customer":
#         flash("Customers may not access this page", "error")
#         return login_manager.unauthorized()
#
#     form = forms.UploadSequencingSampleMappingsForm()
#     from biocomputedm.static.io import sequencing_sample_mappings_template
#     if request.method == "GET":
#         if type == "download":
#             response = excel.make_response_from_book_dict(sequencing_sample_mappings_template.sample_mappings_data,
#                                                           "xls")
#             response.headers["Content-Disposition"] = "attachment; filename=template_for_sequencing_sample_mapping.xls"
#             return response
#         return render_template("input_sequencing_sample_mappings.html", title="Input Sequencing Sample Mappings",
#                                form=form, sid=sid)
#
#     elif type == "upload":
#         if form.validate_on_submit():
#             if form.file_upload.has_file():
#                 filename = request.files["file_upload"].filename
#                 extension = filename.split(".")[1]
#
#                 # Magic to get an csv/excel, transpose it and assign the first value as the dict key
#                 try:
#                     sheet = pyexcel.load_from_memory(extension, request.files["file_upload"].read(),
#                                                      "Sample Mappings Data")
#
#                 except ValueError:
#                     flash("Could not find the Sample Mappings Data sheet in the provided excel file", "error")
#                     return render_template("input_sequencing_sample_mappings.html",
#                                            title="Input Sequencing Sample Mappings", form=form, sid=sid)
#
#                 raw = pyexcel.to_array(sheet)
#                 transposed = list(zip(*raw))
#                 data = dict([(k[0], k[2:]) for k in transposed])
#
#                 passed = sequencing_sample_mappings_template.validate_data_sheet(data)
#                 if not passed:
#                     return render_template("input_sequencing_sample_mappings.html",
#                                            title="Input Sequencing Sample Mappings", form=form, sid=sid)
#
#                 passed = sequencing_sample_mappings_template.build_sample_mappings(sid, data)
#                 if passed:
#                     flash("Document uploaded successfully", "success")
#                     return redirect(url_for("sequencing_submission", sid=sid, type=""))
#
#                 else:
#                     return render_template("input_sequencing_sample_mappings.html",
#                                            title="Input Sequencing Sample Mappings", form=form, sid=sid)
#
#             else:
#                 flash("Missing file information", "error")
#
#         else:
#             flash("Missing file information", "error")
#
#     return render_template("input_sequencing_sample_mappings.html", title="Input Sequencing Sample Mappings", form=form,
#                            sid=sid)
#
#
# @app.route("/run_pipeline/<sample_group>|<pipeline>")
# @login_required("ANY")
# def run_pipeline(sample_group="", pipeline=""):
#     g = utils.get_allowed_sample_group_by_display_key(sample_group)
#     p = utils.get_pipeline_by_display_key(pipeline)
#
#     i = utils.create_pipeline_instance(p)
#     g.current_pipeline = i
#
#     db.session.add(g)
#     db.session.commit()
#
#     # TODO: Start the pipeline!
#
#     return redirect(url_for("sequencing_sample_group", gid=g.display_key))
#
#
# @app.route("/pipeline/<pid>")
# def pipeline(pid=""):
#     p = utils.get_pipeline_by_display_key(pid)
#     return redirect(url_for("empty"))
#
#

#
#
# # TODO Start
# @app.route("/input_flow_cytometry_submission")
# @login_required("ANY")
# def input_flow_cytometry_submission():
#     flash("New flow cytometry submission page is still in development", "warning")
#     return redirect(url_for("empty"))
#
#

