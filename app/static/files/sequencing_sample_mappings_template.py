from app import utils, db
from flask import flash

__author__ = 'jon'

sample_mappings_data = {
    "Sample Mappings Data":
        [
            [
                "Sequencing Samples Group Name",
                "Internal Sample Name",
                "Customer Sample Name",
                "Sequencing Type",
                "Lane Number",
                "Sequencing Concentration",
                "PhiXSpiked",
                "Spike",
                "Spike Ratio",
                "Index 1 Tag Sequence",
                "Index 2 Tag Sequence",
                "Index 1 Tag ID",
                "Index 2 Tag ID",
                "Index 1 Tag Kit ID",
                "Index 2 Tag Kit ID",
                "Adaptor Sequence"
            ],
            [
                "A name for this group of samples, you can have multiple per file or just one. If you leave blank, these samples will be grouped with other samples that are blank for this data submission",
                "Your reference name for this sample",
                "Your customer's reference for this sample - can be edited later!",
                "Type of sequencing. Either rnaseq, chipseq, exomeseq or wgs",
                "The lane for this sample",
                "Default units pm, please specify others, i.e. um or nm if this is not correct",
                "Should be a number!",
                "Should be a number!",
                "Should be a number!",
                "",
                "",
                "Should be a number!",
                "Should be a number!",
                "Should be a number!",
                "Should be a number!",
                ""
            ],
        ],
    "Information":
        [
            ["Core Information"],
            ["If you cannot read"],
            ["All the information"],
            ["below then please"],
            ["widen the columns!"],
            [""],
            ["This file can be used as a template for your data submissions"],
            ["Do NOT change the name of any of the sheets or the column keys"],
            ["Note: Please leave columns 1 and 2 as they are, you can widen them as much as you want!."],
            ["Row 1 contains the key for the information,"],
            ["Row 2 contains a description of how the data should be formatted"]
        ]
}

unessential = [0, 10, 12, 14]

def validate_data_sheet(data):
    d = sample_mappings_data.get("Sample Mappings Data")[0]
    for i in range(0, len(data)):
        if data.get(d[i]) is None:
            flash("Could not identify the column: " + d[i] + " in your submission", "error")
            return False

        if i in unessential:
            continue

        for j in range(0, len(data.get(d[0]))):
            if data.get(d[i])[j] is None or data.get(d[i])[j] == "":
                flash("File was missing data for: " + d[i] + " on line: " + str(
                    j) + " this information is considered essential.", "error")
                return False

    return True


def build_sample_mappings(sid, data):
    d = sample_mappings_data.get("Sample Mappings Data")[0]
    s = utils.get_allowed_submission(sid)
    if s is None:
        flash("Cannot link this sample mappings set to a submission!", "error")
        return False

    for i in range(0, len(data.get(d[0]))):
        if utils.get_allowed_sample_by_internal_name_from_submission(sid, data.get(d[1][i])):
            flash(
                "One of your samples in your mappings file has the same name as another in this submission. This name was: " + data.get(
                    d[1]), "error")
            return False

        t = utils.create_tag(data.get(d[9][i]), data.get(d[11])[i], data.get(d[13])[i], True)
        if t is None:
            return False

        if data.get(d[10])[i] is not None and not data.get(d[10])[i] == "":
            t2 = utils.create_tag(data.get(d[10])[i], data.get(d[12])[i], data.get(d[14])[i], False)
            if t2 is None:
                return False

            sample = utils.create_sequencing_sample(data.get(d[1])[i], data.get(d[2])[i], data.get(d[15])[i], data.get(d[3])[i], t, t2)
        else:
            sample = utils.create_sequencing_sample(data.get(d[1])[i], data.get(d[2])[i], data.get(d[15])[i], data.get(d[3])[i], t)

        s.sample.append(sample)

        sample_group_name = data.get(d[0])[i]
        if sample_group_name is None or sample_group_name == "":
            sample_group_name = "Default"

        g = utils.get_allowed_sample_group_by_name_from_submission(sid, sample_group_name)
        if g is None:
            g = utils.create_sequencing_sample_group(sid, sample_group_name)

        g.sample.append(sample)

        l = utils.get_allowed_lane_by_number_from_submission(sid, data.get(d[4])[i])
        if l is None:
            l = utils.create_lane(sid, data.get(d[4])[i], data.get(d[5])[i], data.get(d[6])[i], data.get(d[7])[i], data.get(d[8])[i])

        sample.lane = l

        db.session.add(sample)
        db.session.add(g)
        db.session.add(l)
        db.session.commit()

    db.session.add(s)
    db.session.commit()

    return True
