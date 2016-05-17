from biocomputedm.admin.models import ReferenceData
from flask.ext.wtf import Form
from flask.ext.wtf.file import FileField, FileRequired
from wtforms import StringField, BooleanField, SelectField, SubmitField
from wtforms.validators import DataRequired


# class RunPipelineForm(Form):
#
#
# class SelectDataForPipelineForm(Form):


class PipelinePropertiesForm(Form):
    execution_field = SelectField("Execution Type:", choices=[("Continuous", "Continuous"), ("Per Module", "Stepwise")])
    options_field = SelectField("Options Choices:", choices=[("Default", "Default"), ("Custom", "Custom")])
    submit_field = SubmitField("Submit")


def build_options_form(options, post_data):
    clazz = OptionsForm
    synthetics = clazz.synthetics
    for synthetic in synthetics:
        clazz = clazz.remove_field(synthetic)
    synthetics.clear()

    # Loop through options and append a form field based on type
    for option in options:
        if option.user_interaction_type == "file":
            clazz = clazz.append_field(option.display_key + "_upload",
                                       FileField(option.display_name,
                                                 validators=[FileRequired()],
                                                 description=option.description))
            synthetics.append(option.display_key + "_upload")

            # Optional template provider - should be default though!
            if option.default_value is not None and option.default_value != "":
                clazz = clazz.append_field(option.display_key + "_template",
                                           SubmitField("Download Template:",
                                                       description="Provides a template for the required file upload"))
                synthetics.append(option.display_key + "_template")

        elif option.user_interaction_type == "string":
            clazz = clazz.append_field(option.display_key,
                                       StringField(option.display_name,
                                                   default=option.default_value,
                                                   validators=[DataRequired()],
                                                   description=option.description))
            synthetics.append(option.display_key)

        elif option.user_interaction_type == "boolean":
            clazz = clazz.append_field(option.display_key,
                                       BooleanField(option.display_name,
                                                    description=option.description))
            synthetics.append(option.display_key)

        elif option.user_interaction_type == "reference":
            choices = []
            choice_options = ReferenceData.query.filter_by(current=True).all()
            for choice_option in choice_options:
                choices.append((choice_option.display_key, choice_option.name + " (" + choice_option.version + "): " + choice_option.description))

            if len(choices) == 0:
                choices.append(("NO", "No libraries found..."))

            clazz = clazz.append_field(option.display_key,
                                       SelectField(option.display_name,
                                                   choices=choices,
                                                   validators=[DataRequired()],
                                                   description=option.description))

            synthetics.append(option.display_key)

        elif option.user_interaction_type == "enum":
            choose = option.default_value.split(",")
            choices = []
            for choice in choose:
                choices.append((choice, choice))

            clazz = clazz.append_field(option.display_key,
                                       SelectField(option.display_name,
                                                   default=choices[0],
                                                   choices=choices,
                                                   validators=[DataRequired()],
                                                   description=option.description))
            synthetics.append(option.display_key)

        else:
            continue

    clazz.synthetics = synthetics
    form = clazz()
    # if post_data is not None:
    #     form = clazz(post_data)
    # else:
    #     form = clazz()

    return form


class OptionsForm(Form):
    synthetics = []
    submit_field = SubmitField("Submit")

    @classmethod
    def append_field(cls, name, field):
        setattr(cls, name, field)
        return cls

    @classmethod
    def remove_field(cls, name):
        delattr(cls, name)
        return cls
