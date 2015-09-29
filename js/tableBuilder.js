/**
 * Created by jon on 28/09/15.
 */

function buildProjectsTable(projects) {
    var json = JSON.parse(projects);
    var table = document.getElementById('projectstable');
    var body = table.getElementsByTagName('tbody')[0];
    for (var i = 0; i < json.length; i++) {
        var newRow = body.insertRow(body.rows.length);

        for (var j = 0; j < 7; j++) {
            var newCell = newRow.insertCell(j);

            var newNode;
            switch (j) {
                case 0:
                    newNode = document.createTextNode(json[i].masterProjectID);
                    break;

                case 1:
                    newNode = document.createElement("a");
                    newNode.setAttribute("href", "/cgi-bin/getProjectFromDB?projName=" + json[i].projectName + "&projID=" + json[i].masterProjectID);
                    var childText = document.createTextNode(json[i].projectName);
                    newNode.appendChild(childText);
                    break;

                case 2:
                    newNode = document.createTextNode(json[i].projectLead);
                    break;

                case 3:
                    newNode = document.createTextNode("Customer");
                    break;

                case 4:
                    newNode = document.createTextNode(json[i].openDate);
                    break;

                case 5:
                    newNode = document.createTextNode(json[i].lastUpdate);
                    break;

                case 6:
                    newNode = document.createTextNode(json[i].status);
                    break;
            }

            newCell.appendChild(newNode);
        }
    }

    //for pro in range(0, len(projects)):
    //    # print "<p>",projects[pro],"</p>"
    //        if len(projects[pro]) > 7:
    //        print \
    //        '<tr><td>' + str(projects[pro][0]) + '</td><td><a href="/cgi-bin/getProjectFromDB?projName=' + \
    //        projects[pro][1] + '&projID=' + str(projects[pro][0]) + '"> ' + projects[pro][1] + '</a></td><td>' + \
    //        projects[pro][2] + '</td><td>Customer</td><td>' + str(projects[pro][5]) + '</td><td>' + str(
    //            projects[pro][6]) + '</td><td>' + projects[pro][3] + '</td></tr>'
    //    else:
    //        print \
    //        '<tr><td>' + str(projects[pro][0]) + '</td><td><a href="/cgi-bin/getProjectFromDB?projName=' + \
    //        projects[pro][1] + '&projID=' + str(projects[pro][0]) + '"> ' + projects[pro][1] + '</a></td><td>' + \
    //        projects[pro][2] + '</td><td>Customer</td><td>' + str(projects[pro][5]) + '</td><td>' + str(
    //            projects[pro][6]) + '</td><td>' + projects[pro][3] + '</td></tr>'
}

function buildSequenceRunsTable(sequencingRuns) {
    var json = JSON.parse(sequencingRuns);
    var table = document.getElementById('sequencerunstable');
    var body = table.getElementsByTagName('tbody')[0];
    for (var i = 0; i < json.length; i++) {
        var newRow = body.insertRow(body.rows.length);

        for (var j = 0; j < 5; j++) {
            var newCell = newRow.insertCell(j);

            var newNode;
            switch (j) {
                case 0:
                    newNode = document.createElement("a");
                    newNode.setAttribute("href", "/cgi-bin/getSeqRunFromDB?projID=" + json[i].seqRunID);
                    var childText = document.createTextNode(json[i].seqRunName);
                    newNode.appendChild(childText);
                    break;

                case 1:
                    newNode = document.createTextNode(json[i].startDate);
                    break;

                case 2:
                    newNode = document.createTextNode(json[i].completionDate);
                    break;

                case 3:
                    newNode = document.createTextNode(json[i].seqProjectName);
                    break;

                case 4:
                    if (json[i].masterProjectID != null && json[i].projectName != null) {
                        newNode = document.createElement("a");
                        newNode.setAttribute("href", "/cgi-bin/getProjectFromDB?projName=" + json[i].projectName + "&projID=" + json[1].masterProjectID);
                        var childText = document.createTextNode(json[i].projectName);
                        newNode.appendChild(childText);
                        break;
                    }

                    else {
                        newNode = document.createTextNode(json[i].projectName);
                    }
            }

            newCell.appendChild(newNode);
        }
    }

    //for exp in range(0, len(experiments)):
    //for pro in range(0, len(experiments[exp])):
    //print \
    //        '<tr><td><a href="/cgi-bin/getSeqRunFromDB?projID=' + str(experiments[exp][pro][0]) + '"> ' + str(
    //            experiments[exp][pro][3]) + ' </a></td><td>' + str(experiments[exp][pro][1]) + '</td><td>' + str(
    //            experiments[exp][pro][2]) + '</td><td>' + str(
    //            experiments[exp][pro][4]) + '</td><td><a href="/cgi-bin/getProjectFromDB?projName=' + str(
    //            experiments[exp][pro][6]) + '&projID=' + str(experiments[exp][pro][5]) + '">' + str(
    //            experiments[exp][pro][6]) + '</a></td></tr>'
}
