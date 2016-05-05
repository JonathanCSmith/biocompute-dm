from biocomputedm import utils
from biocomputedm.admin import forms
from biocomputedm.admin import models
from biocomputedm.admin.models import User, ReferenceData, Person, Customer, UserGroup, CustomerGroup, Group
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
            flash("Successfully logged in.", "success")
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
            UserGroup.create(
                    group_name=str(form.group_name.data),
                    admin_name=str(form.admin_login.data),
                    admin_password=str(form.admin_password.data),
                    admin_email=str(form.admin_email.data)
            )
            return redirect(url_for("admin.show_groups"))


@admin.route("/show_users")
@admin.route("/show_users/<int:page>")
@login_required("Site Admin", "Group Admin")
def show_users(page=1):
    if current_user.get_role() == "Site Admin":
        u = Person.query.paginate(page=page, per_page=20)

    elif current_user.type == "Customer":
        u = current_user.group.members.paginate(page=page, per_page=20)

    else:
        u = current_user.group.members.paginate(page=page, per_page=20)

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
    g = current_user.group.customer_groups
    if g is not None:
        g = g.paginate(page=page, per_page=20)
    return render_template("groups.html", title="Consignors", page=page, obs=g)


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
            customer_group = CustomerGroup.create(
                group_name=str(form.group_name.data),
                admin_name=str(form.admin_login.data),
                admin_password=str(form.admin_password.data),
                admin_email=str(form.admin_email.data),
                parent_group=current_user.group
            )

            return redirect(url_for("admin.show_customers"))


@admin.route("/add_customer", methods=["GET", "POST"])
@login_required("Group Admin")
def add_customer():
    form = forms.CreatePerson()
    if request.method == "GET":
        return render_template("add_person.html", title="Add Consignor", type="User", form=form)

    else:
        if not form.validate():
            utils.flash_errors(form)
            return render_template("add_person.html", title="Add Consignor", type="User", form=form)

        else:
            Customer.create(
                username=str(form.login_name.data),
                password=str(form.login_password.data),
                email=str(form.login_email.data),
                group=current_user.group
            )
            return redirect(url_for("admin.show_users"))


@admin.route("/link_to_customer/<oid>|<origin>", methods=["GET", "POST"])
@login_required("ANY")
def link_to_customer(oid="", origin=""):
    if oid == "" or origin != "project":
        flash("Could not locate the provided information", "error")
        return redirect(url_for("index"))

    project = current_user.group.projects.filter_by(display_key=oid).first()
    if project is None:
        flash("Could not identify the provided project.", "error")
        return redirect(url_for("index"))

    from biocomputedm.admin import forms
    form = forms.UpdateCustomerLinkForm()
    if request.method == "POST":
        # Check to see that at least 1 upload was selected - empty submissions have no use
        ids = request.form.getlist("do_select")
        if ids is not None and len(ids) != 0:
            for key in ids:
                group = CustomerGroup.query.filter_by(display_key=key).first()
                if group is not None:
                    group.projects.append(project)
                    group.save()

                    for sample in project.samples:
                        group.samples.append(sample)

                    group.save()

                    for pipeline_output in project.pipeline_outputs:
                        group.data_groups.append(pipeline_output)

                    group.save()

            flash("Consignor was successfully updated", "success")
            return redirect(url_for("index"))

        else:
            flash("No consignors were selected.", "warning")
            return redirect(url_for("index"))

    consignors = CustomerGroup.query.filter_by(parent_id=current_user.group.id).all()
    potential_consignors = []
    for consignor in consignors:
        skip = False
        projects = consignor.projects
        for project in projects:
            if project.id == project.id:
                skip = True
                break

        if skip:
            continue

        else:
            potential_consignors.append(consignor)

    return render_template("select_customer.html", title="Link to Consignor", customers=consignors, obj=project, origin=origin, form=form)


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
