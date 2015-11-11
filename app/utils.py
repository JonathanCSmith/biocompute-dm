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


def get_allowed_sample_by_internal_name_from_submission(sid, name):
    return get_allowed_sequencing_samples_from_submission_query(sid).filter_by(internal_sample_name=name).first()


def get_allowed_sequencing_samples_from_submission_query(sid):
    return SequencingSample.query.filter_by(submission_id=sid)


def create_lane(sid, number, concentration, phi, spike, ratio):
    if current_user.type == "Customer":
        return None

    l = get_allowed_lane_by_number_from_submission(sid, number)
    if l is not None:
        return l

    l = Lane()
    l.number = number
    l.set_sequencing_concentration(concentration)
    l.phi_x_spiked = phi
    l.spike = spike
    l.spike_ratio = ratio

    s = get_allowed_submission(sid)
    s.lane.append(l)

    db.session.add(l)
    db.session.commit()

    return l


def get_allowed_lane_by_number_from_submission(sid, number):
    return Lane.query.filter_by(submission_id=sid).filter_by(number=number).first()


def create_sequencing_sample_group(sid, name):
    if current_user.type == "Customer":
        return None

    s = get_allowed_submission(sid)
    g = get_allowed_sample_group_by_name_from_submission(sid, name)
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


def get_allowed_sample_group_by_name_from_submission(sid, name):
    return SampleGroup.query.filter_by(submission_id=sid, name=name).first()


def get_allowed_sample_group(gid):
    if current_user.is_authenticated:
        return SampleGroup.query.filter_by(group_id=current_user.group_id, id=gid).first()


def get_allowed_sample_groups_query():
    if current_user.is_authenticated:
        return SampleGroup.query.filter_by(group_id=current_user.group_id)


def create_sequencing_submission(name, flowid, start_date, end_date, lead, location, cyclest1, cyclest2, cyclesr1, cyclesr2, pe):
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


def get_allowed_submission(sid):
    return get_allowed_submissions_query().filter_by(id=sid).first()


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


def link_sample_group_to_investigation(iid, gid):
    i = get_allowed_investigation(iid)
    p = get_allowed_sample_group(gid)

    if i is None or p is None:
        return False

    i.sample_group.append(p)
    db.session.add(i)
    db.session.commit()
    return True


def get_allowed_investigation(iid):
    return get_allowed_investigation_query(iid).first()


def get_allowed_investigation_query(iid):
    i = get_allowed_investigations_query()
    return i.filter_by(id=iid)


def get_allowed_investigations():
    return get_allowed_investigations_query().all()


def get_allowed_investigations_query():
    if current_user.is_authenticated:
        return Investigation.query.filter_by(group_id=current_user.group_id)

    return None


def create_document(iid, filename, description, filepath):
    i = get_allowed_investigation(iid)
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


def remove_document(did):
    doc = get_allowed_document(did)
    i = get_allowed_investigation(doc.investigation_id)
    i.document.remove(doc)
    db.session.delete(doc)
    db.session.commit()


def get_allowed_document(did):
    return get_allowed_documents_query().filter_by(id=did).first()


def get_allowed_documents_query():
    if current_user.is_authenticated:
        return Document.query.filter_by(group_id=current_user.group_id)

    return None