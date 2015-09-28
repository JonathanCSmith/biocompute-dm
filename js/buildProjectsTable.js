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
                    newNode = document.createTextNode(json[i][0]);
                    break;

                case 1:
                    newNode = document.createElement("a");
                    newNode.setAttribute("href", "/cgi-bin/getProjectFromDB?projName=" + json[i][1] + "&projID=" + json[i][0]);
                    var childText = document.createTextNode(json[i][1]);
                    newNode.appendChild(childText);
                    break;

                case 2:
                    newNode = document.createTextNode(json[i][2]);
                    break;

                case 3:
                    newNode = document.createTextNode("Customer");
                    break;

                case 4:
                    newNode = document.createTextNode(json[i][5]);
                    break;

                case 5:
                    newNode = document.createTextNode(json[i][6]);
                    break;

                case 6:
                    newNode = document.createTextNode(json[i][3]);
                    break;
            }

            newCell.appendChild(newNode);
        }
    }
}

/**
for pro in range(0, len(projects)):
# print "<p>",projects[pro],"</p>"
if len(projects[pro]) > 7:
print \
    '<tr><td>' + str(projects[pro][0]) + '</td><td><a href="/cgi-bin/getProjectFromDB?projName=' + \
    projects[pro][1] + '&projID=' + str(projects[pro][0]) + '"> ' + projects[pro][1] + '</a></td><td>' + \
    projects[pro][2] + '</td><td>Customer</td><td>' + str(projects[pro][5]) + '</td><td>' + str(
        projects[pro][6]) + '</td><td>' + projects[pro][3] + '</td></tr>'
else:
print \
    '<tr><td>' + str(projects[pro][0]) + '</td><td><a href="/cgi-bin/getProjectFromDB?projName=' + \
    projects[pro][1] + '&projID=' + str(projects[pro][0]) + '"> ' + projects[pro][1] + '</a></td><td>' + \
    projects[pro][2] + '</td><td>Customer</td><td>' + str(projects[pro][5]) + '</td><td>' + str(
        projects[pro][6]) + '</td><td>' + projects[pro][3] + '</td></tr>'
 **/