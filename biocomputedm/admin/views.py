from biocomputedm import utils
from biocomputedm.admin import forms
from biocomputedm.admin.models import User, Group, create_group, Customer
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
            user = User.query.filter_by(username=str(form.username.data)).first()
            flash("Successfully logged in!", "success")
            login_user(user)
            return redirect(url_for("index"))


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
@login_required("Site Admin", "Group Admin")
def show_customers(page=1):
    if current_user.get_role() == "Site Admin":
        c = Customer.query.paginate(page=page, per_page=20)

    else:
        g = current_user.group
        c = Customer.query.filter_by(group_id=g.id).paginate(page=page, per_page=20)

    return render_template("people.html", title="Customers", page=page, obs=c)


@admin.route("/add_customer", methods=["GET", "POST"])
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
            Customer.create(
                    username=str(form.login_name.data),
                    password=str(form.login_password.data),
                    email=str(form.login_password),
                    group=current_user.group
            )
            return redirect(url_for("admin.show_customers"))
