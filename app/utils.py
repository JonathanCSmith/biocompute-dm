import config
from app.models import *
from flask import flash
from flask.ext.login import current_user

__author__ = 'jon'


# TODO: We could stop duplication here if we are smart
def create_tag(tag_sequence, tag_id, tag_kit_id, isfirst):
    if tag_sequence is None and (tag_id is None or tag_kit_id is None):
        flash("There was an error processing your tag information", "error")

    elif tag_sequence is None:
        tag_sequence = ""

    elif tag_id is None and tag_kit_id is None:
        tag_id = -1
        tag_kit_id = -1

    elif tag_id is None:
        tag_id = -1

    elif tag_kit_id is None:
        tag_kit_id = -1

    t = Tag()
    t.tag_id = tag_id
    t.tag_library = tag_kit_id
    t.tag_sequence = tag_sequence
    t.is_first_index = isfirst

    db.session.add(t)
    db.session.commit()

    return t


def create_sequencing_sample(internal_name, customer_name, adaptor_sequence, sample_type, tag1, tag2=None):
    if current_user.type == "Customer":
        return

    s = SequencingSample()
    s.internal_sample_name = internal_name
    s.customer_sample_name = customer_name
    s.adaptor_sequence = adaptor_sequence
    s.index_tag.append(tag1)
    s.sample_type = sample_type
    if tag2 is not None:
        s.index_tag.append(tag2)

    current_user.group.sample.append(s)
    current_user.sample.append(s)

    db.session.add(s)
    db.session.add(current_user.group)
    db.session.add(current_user)
    db.session.commit()

    return s


def get_allowed_sample_by_internal_name_from_submission_display_key(sample_display_key, name):
    return get_allowed_sequencing_samples_from_submission_display_key_query(sample_display_key).filter_by(
        internal_sample_name=name).first()


def get_allowed_sequencing_samples_from_submission_display_key_query(display_key):
    sid = Submission.query.filter_by(display_key=display_key).first().id
    return SequencingSample.query.filter_by(submission_id=sid)


def create_lane(sid, number, concentration, phi, spike, ratio):
    if current_user.type == "Customer":
        return None

    l = get_allowed_lane_by_number_from_submission_display_key(sid, number)
    if l is not None:
        return l

    l = Lane()
    l.number = number
    l.set_sequencing_concentration(concentration)
    l.phi_x_spiked = phi
    l.spike = spike
    l.spike_ratio = ratio

    s = get_allowed_submission_by_display_key(sid)
    s.lane.append(l)

    db.session.add(l)
    db.session.commit()

    return l


def get_allowed_lane_by_number_from_submission_display_key(display_key, number):
    sid = Submission.query.filter_by(display_key=display_key).first().id
    return Lane.query.filter_by(submission_id=sid).filter_by(number=number).first()


def create_sequencing_sample_group(sid, name):
    if current_user.type == "Customer":
        return None

    s = get_allowed_submission_by_display_key(sid)
    g = get_allowed_sample_group_by_name_from_submission_display_key(sid, name)
    if g is not None:
        return g

    g = SequencingSampleGroup()
    g.name = name

    current_user.group.sample_group.append(g)
    current_user.sample_group.append(g)
    s.sample_group.append(g)

    db.session.add(s)
    db.session.add(g)
    db.session.add(current_user.group)
    db.session.add(current_user)
    db.session.commit()

    return g


def get_allowed_sample_group_by_name_from_submission_display_key(display_key, name):
    sid = Submission.query.filter_by(display_key=display_key).first().id
    return SampleGroup.query.filter_by(submission_id=sid, name=name).first()


def get_allowed_sample_group_by_display_key(display_key):
    if current_user.is_authenticated:
        return SampleGroup.query.filter_by(group_id=current_user.group_id, display_key=display_key).first()


def get_allowed_sample_groups_query():
    if current_user.is_authenticated:
        return SampleGroup.query.filter_by(group_id=current_user.group_id)


def create_sequencing_submission(name, flowid, start_date, end_date, lead, location, cyclest1, cyclest2, cyclesr1,
                                 cyclesr2, pe):
    if current_user.type == "Customer":
        return

    from datetime import date
    if not isinstance(start_date, date) or not isinstance(end_date, date):
        flash("Incorrect format for date time in submission", "error")
        return None

    if not isinstance(cyclest1, int) and not isinstance(cyclest1, float):
        flash("Incorrect format for the number of tag cycles", "error")
        return None

    if not isinstance(cyclest2, int) and not isinstance(cyclest2, float):
        flash("Incorrect format for the number of tag 2 cycles", "error")
        return None

    if not isinstance(cyclest1, int) and not isinstance(cyclest1, float):
        flash("Incorrect format for the number of reads for 1", "error")
        return None

    if not isinstance(cyclest1, int) and not isinstance(cyclest1, float):
        flash("Incorrect format for the number of reads for 2", "error")
        return None

    if not pe == "Yes" and not pe == "No":
        flash("Incorrect format for the paired end information", "error")
        return None

    s = SequencingSubmission()
    s.name = name
    s.flow_cell_id = flowid
    s.start_date = start_date
    s.completion_date = end_date
    s.leader = lead
    s.data_location = location
    s.index_tag_cycles = cyclest1
    s.index_tag_cycles_2 = cyclest2
    s.read_cycles = cyclesr1
    s.read_cycles_2 = cyclesr2
    s.paired_end = pe

    current_user.group.submission.append(s)
    current_user.submission.append(s)

    db.session.add(current_user)
    db.session.add(current_user.group)
    db.session.add(s)
    db.session.commit()
    return s


def get_allowed_submission_by_display_key(display_key):
    return get_allowed_submissions_query().filter_by(display_key=display_key).first()


def get_allowed_submissions_query():
    if current_user.is_authenticated:
        return Submission.query.filter_by(group_id=current_user.group_id)


def create_investigation(name, lead):
    i = Investigation(name, lead)

    current_user.investigation.append(i)
    current_user.group.investigation.append(i)

    db.session.add(i)
    db.session.add(current_user)
    db.session.commit()
    return i


def link_sample_group_to_investigation(investigation_display_key, group_display_key):
    i = get_allowed_investigation_by_display_key(investigation_display_key)
    p = get_allowed_sample_group_by_display_key(group_display_key)

    if i is None or p is None:
        return False

    i.sample_group.append(p)
    db.session.add(i)
    db.session.commit()
    return True


def get_allowed_investigation_by_display_key(display_key):
    return get_allowed_investigation_query_by_display_key(display_key).first()


def get_allowed_investigation_by_id(id):
    return get_allowed_investigation_query_by_id(id).first()


def get_allowed_investigation_query_by_display_key(display_key):
    i = get_allowed_investigations_query()
    return i.filter_by(display_key=display_key)


def get_allowed_investigation_query_by_id(id):
    i = get_allowed_investigations_query()
    return i.filter_by(id=id)


def get_allowed_investigations():
    return get_allowed_investigations_query().all()


def get_allowed_investigations_query():
    if current_user.is_authenticated:
        return Investigation.query.filter_by(group_id=current_user.group_id)

    return None


def create_document(investigation_display_key, filename, description, filepath):
    i = get_allowed_investigation_by_display_key(investigation_display_key)
    if i is None:
        return

    doc = Document()
    doc.name = filename
    doc.description = description
    doc.location = filepath

    current_user.document.append(doc)
    current_user.group.document.append(doc)
    i.document.append(doc)

    db.session.add(i)
    db.session.add(doc)
    db.session.add(current_user)
    db.session.commit()
    return doc


def remove_document(document_display_key):
    doc = get_allowed_document_by_display_key(document_display_key)
    i = get_allowed_investigation_by_display_key(doc.investigation_id)
    i.document.remove(doc)
    db.session.delete(doc)
    db.session.commit()


def get_allowed_document_by_display_key(document_display_key):
    return get_allowed_documents_query().filter_by(display_key=document_display_key).first()


def get_allowed_documents_query():
    if current_user.is_authenticated:
        return Document.query.filter_by(group_id=current_user.group_id)

    return None


def get_allowed_pipeline(name, description, author, version):
    if current_user.is_authenticated and current_user.role == "Site Admin":
        return Pipeline.query.filter_by(name=name, description=description, author=author, version=version).first()

    return None


def get_pipelines():
    return Pipeline.query.all()


def get_pipeline_by_display_key(display_key):
    return Pipeline.query.filter_by(display_key=display_key).first()


def create_pipeline(name, description, author, version):
    pipeline = Pipeline()
    pipeline.name = name
    pipeline.description = description
    pipeline.author = author
    pipeline.version = version

    db.session.add(pipeline)
    db.session.commit()
    return pipeline


def create_module(name, description, executor, order_index, pipeline):
    module = PipelineModule()
    module.name = name
    module.description = description

    # Full path to executor - we have to use the config value here because it may be located at a different place
    # on the remote
    p = os.path.join(config.REMOTE_PIPELINES_PATH, pipeline.name)
    p = os.path.join(p, executor)
    module.executor = p

    module.execution_index = order_index
    pipeline.module.append(module)

    db.session.add(module)
    db.session.add(pipeline)
    db.session.commit()
    return module


def create_option(display_name, parameter_name, default_value, user_interaction_type, module):
    option = PipelineModuleOption()
    option.display_name = display_name
    option.paramater_name = parameter_name
    option.user_interaction_type = user_interaction_type
    option.default_value = default_value
    module.module_option.append(option)

    db.session.add(option)
    db.session.add(module)
    db.session.commit()
    return option


def create_pipeline_instance(p):
    pipeline = PipelineInstance()
    p.instance.append(pipeline)

    db.session.add(pipeline)
    db.session.add(p)
    db.session.commit()

    return pipeline


# Function to display errors
def flash_errors(form):
    for field, errors in form.errors.items():
        for error in errors:
            flash(u"Error in the %s field - %s" % (
                getattr(form, field).label.text,
                error
            ), "error")


def refresh_pipelines():
    path = os.path.join(os.path.dirname(__file__), "..")
    path = os.path.join(path, "link")
    if not os.path.exists(path):
        flash("Could not locate link directory. Please ensure this folder exists.", "error")
        return False

    path = os.path.join(path, "pipelines")
    if not os.path.exists(path):
        flash("Could not locate pipelines directory. Please ensure there is a valid symlink for the folder in the links folder!", "error")
        return False

    directories = os.listdir(path)
    found = False
    for directory in directories:
        directory_path = os.path.join(path, directory)
        if not os.path.isdir(directory_path):
            continue

        file = os.path.join(directory_path, directory + ".json")
        if not os.path.isfile(file):
            continue

        from app.static.io import pipeline_mappings_template as template_helper
        if not template_helper.validate(file):
            continue

        found |= template_helper.build(file)

    if found:
        flash("Successfully loaded all pipelines from the pipeline directory.", "success")
    else:
        flash("No pipelines were loaded as none were found or they were already present.", "warning")

    return True