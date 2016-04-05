from biocomputedm import utils
from biocomputedm.admin import forms
from biocomputedm.admin import models
from biocomputedm.admin.models import User, Group, create_group, ReferenceData, create_customer_group, Person
from biocomputedm.decorators import login_required
from flask import Blueprint, request, render_template, redirect, url_for, flash, g
from flask.ext.login import login_user, logout_user, current_user

admin = Blueprint("admin", __name__, static_folder="static", template_folder="templates")


@admin.route("/login", methods=["GET", "POST"])
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
            user = Person.query.filter_by(username=str(form.username.data)).first()
            flash("Successfully logged in!", "success")
            login_user(user)

            next = request.args.get('next')
            return redirect(next or url_for('index'))


@admin.route("/logout")
@login_required("ANY")
def logout():
    logout_user()
    return redirect(url_for("index"))


@admin.route("/administrate")
@login_required("Site Admin", "Group Admin")
def administrate():
    return render_template("administrate_welcome.html", title="Administrate Panel")


@admin.route("/show_groups")
@admin.route("/show_groups/<int:page>")
@login_required("Site Admin")
def show_groups(page=1):
    g = Group.query.paginate(page=page, per_page=20)
    return render_template("groups.html", title="Groups", page=page, obs=g)


@admin.route("/add_group", methods=["GET", "POST"])
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
            create_group(
                    str(form.group_name.data),
                    str(form.admin_login.data),
                    str(form.admin_password.data),
                    str(form.admin_email.data)
            )
            return redirect(url_for("admin.show_groups"))


@admin.route("/show_users")
@admin.route("/show_users/<int:page>")
@login_required("Site Admin", "Group Admin")
def show_users(page=1):
    if current_user.get_role() == "Site Admin":
        u = User.query.paginate(page=page, per_page=20)

    else:
        g = current_user.group
        u = User.query.filter_by(group_id=g.id).paginate(page=page, per_page=20)

    return render_template("people.html", title="Users", page=page, obs=u)


@admin.route("/add_user", methods=["GET", "POST"])
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
            User.create(
                    username=str(form.login_name.data),
                    password=str(form.login_password.data),
                    email=str(form.login_email.data),
                    group=current_user.group
            )
            return redirect(url_for("admin.show_users"))


@admin.route("/show_customers")
@admin.route("/show_customers/<int:page>")
@login_required("Group Admin")
def show_customers(page=1):
    g = Group.query.filter_by(parent_id=current_user.group.id)
    if g is not None:
        g = g.paginate(page=page, per_page=20)
    return render_template("groups.html", title="Consigner", page=page, obs=g)


@admin.route("/add_customer_group", methods=["GET", "POST"])
@login_required("Group Admin")
def add_customer_group():
    form = forms.CreateGroupForm()
    if request.method == "GET":
        return render_template("add_group.html", title="Consignor Creation", form=form)

    else:
        if not form.validate():
            utils.flash_errors(form)
            return render_template("add_group.html", title="Consignor Creation", form=form)

        else:
            customer_group = create_customer_group(
                    str(form.group_name.data),
                    str(form.admin_login.data),
                    str(form.admin_password.data),
                    str(form.admin_email.data)
            )

            customer_group.parent_id = current_user.group.id
            customer_group.save()

            return redirect(url_for("admin.show_customers"))


@admin.route("/link_to_customer/<oid>|<origin>", methods=["GET", "POST"])
@login_required("ANY")
def link_to_customer(oid="", origin=""):
    if oid == "" or origin == "":
        flash("Could not locate the provided sample group", "error")
        return redirect(url_for("index"))

    if origin == "sample_group":
        obj = current_user.group.sample_groups.filter_by(display_key=oid).first()
    elif origin == "sample":
        obj = current_user.group.samples.filter_by(display_key=oid).first()
    elif origin == "project":
        obj = current_user.group.projects.filter_by(display_key=oid).first()
    else:
        flash("Could not identify the object to link", "error")
        return redirect(url_for("index"))

    if obj is None:
        flash("Could not identify the object to link", "error")
        return redirect(url_for("index"))

    from biocomputedm.admin import forms
    form = forms.UpdateCustomerLinkForm()
    if request.method == "POST":
        # Check to see that at least 1 upload was selected - empty submissions have no use
        ids = request.form.getlist("do_select")
        if ids is not None and len(ids) != 0:
            for key in ids:
                group = Group.query.filter_by(display_key=key).first()
                if group is not None:
                    if origin == "sample_group":
                        group.sample_groups.append(obj)
                        group.save()

                    elif origin == "sample":
                        group.samples.append(obj)
                        group.save()

                    else:
                        group.projects.append(obj)
                        group.save()

            flash("Consignor was successfully updated", "success")

        else:
            flash("No consignors were selected.", "warning")

    consignors = Group.query.filter_by(parent_id=current_user.group.id).all()
    potential_consignors = []
    for consignor in consignors:
        skip = False
        if origin == "sample_group":
            sample_groups = consignor.sample_groups
            for sample_group in sample_groups:
                if sample_group.id == obj.id:
                    skip = True
                    break

            if skip:
                continue

            else:
                potential_consignors.append(consignor)

        elif origin == "sample":
            samples = consignor.samples
            for sample in samples:
                if sample.id == obj.id:
                    skip = True
                    break

            if skip:
                continue

            else:
                potential_consignors.append(consignor)

        else:
            projects = consignor.projects
            for project in projects:
                if project.id == obj.id:
                    skip = True
                    break

            if skip:
                continue

            else:
                potential_consignors.append(consignor)

    return render_template("select_customer.html", title="Link to Consignor", customers=consignors, obj=obj, origin=origin, form=form)


@admin.route("/refresh_reference_data")
@login_required("Site Admin")
def refresh_reference_data():
    found = models.refresh_reference_data_library()
    if found:
        flash("Successfully loaded all reference data libraries", "success")
    else:
        flash("No reference data libraries were loaded as no new members were identified.", "warning")

    return redirect(url_for("admin.administrate"))


@admin.route("/display_reference_data")
@admin.route("/display_reference_data/<int:page>")
@login_required("ANY")
def display_reference_data(page=1):
    items = ReferenceData.query.paginate(page=page, per_page=20)
    return render_template("reference_libraries.html", title="Reference Data", page=page, obs=items)
