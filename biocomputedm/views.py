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
