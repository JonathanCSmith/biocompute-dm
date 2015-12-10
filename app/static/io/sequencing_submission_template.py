from app import utils

__author__ = 'jon'

submission_data = {
    "Submission Data":
        [
            ["Sequencing Submission Name", "Text used to identify this data submission", ""],
            ["Flow Cell ID", "Unique identifier used to demultiplex and identify information", ""],
            ["Start Date", "When this submission was started. Format according to: YYYY-MM-DD", ""],
            ["Completion Date", "When this run was finished. Format according to: YYYY-MM-DD", ""],
            ["Genomics Lead", "Individual in charge of this submission", ""],
            ["Data Location", "Where your data is located on the gpfs. This will change in future!", ""],
            ["Number of Cycles for Index Tag 1", "Number only!", ""],
            ["Number of Cycles for Index Tag 2", "Number only!"],
            ["Number of Cycles for Read 1", "Number only!", ""],
            ["Number of Cycles for Read 2", "Number only!", ""],
            ["Paired End", "Format: Yes or No", ""]
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
            ["Do NOT change the name of any of the sheets or the row keys"],
            ["Note: Please leave columns 1 and 2 as they are, you can widen them as much as you want!."],
            ["Column 1 contains the key for the information,"],
            ["Column 2 contains a description of how the data should be formatted"]
        ]
}


def validate_data_sheet(data):
    d = submission_data.get("Submission Data")
    for i in range(0, len(data)):
        if not d[i][0] in data:
            return d[i][0]

        if data.get(d[i][0]) == "":
            return d[i][0]

    return None


def build_submission_entry(data):
    d = submission_data.get("Submission Data")
    return utils.create_sequencing_submission(
        data.get(d[0][0]),
        data.get(d[1][0]),
        data.get(d[2][0]),
        data.get(d[3][0]),
        data.get(d[4][0]),
        data.get(d[5][0]),
        data.get(d[6][0]),
        data.get(d[7][0]),
        data.get(d[8][0]),
        data.get(d[9][0]),
        data.get(d[10][0]))
