from biocomputedm.decorators import login_required
from biocomputedm.pipelines.views import pipelines
from flask import Blueprint, render_template, g
from flask import abort

content = Blueprint("content", __name__, static_folder="static", template_folder="templates")


@content.route("/welcome")
def welcome():
    return render_template("welcome.html", title="Welcome", user=g.user)


@content.route("/activity")
@login_required("ANY")
def activity():
    return render_template("activity.html", title="Task Panel")


@content.route("/about")
def about():
    return render_template("welcome.html", title="About Biocompute-DM", user=g.user)


@content.route("/data_processing")
def data_processing():
    return render_template("data_processing.html", title="Data Processing", user=g.user)


@content.route("/data_management")
def data_management():
    return render_template("data_management.html", title="Data Management", user=g.user)


@content.route("/data_monitoring")
def data_monitoring():
    return render_template("data_monitoring.html", title="Data Monitoring", user=g.user)


@content.route("/data_analysis")
def data_analysis():
    return render_template("data_analysis.html", title="Data Assessment", user=g.user)


@content.route("/terms_and_conditions")
def terms_and_conditions():
    return render_template("terms_and_conditions.html", title="Terms and Conditions")


@content.route("/message/<type>|<oid>", methods=["POST"])
def message(type="", oid=""):
    if type == "manage":
        from biocomputedm.manage.views import message
        return message(oid=oid)

    elif type == "pipelines":
        from biocomputedm.pipelines.views import message
        return pipelines.message(oid=oid)

    else:
        return abort(404)
