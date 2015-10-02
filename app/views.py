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


@app.route("/projects")
@login_required("ANY")
def projects():
    flash("Projects page is still in development")
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
    return redirect(url_for("index"))

# @login_required to secure
