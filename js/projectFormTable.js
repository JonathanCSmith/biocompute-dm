function experimentTable(exptArr) {

    alert("In expt table");
    function addSample(lastSample) {

        var sampleRow = document.getElementById("sampleRow_" + lastSample);
        var laneCell1 = document.getElementById("laneCell_" + lastSample);
        var sampRows = laneCell1.getElementsByTagName('input');

        //alert(sampRows[sampRows.length-1].name);

        //alert(sampleRow.id);
        //alert(laneCell1.id);
        var newRow = sampleRow.cloneNode(true);
        var newRowInputs = newRow.getElementsByTagName('input');
        //alert(newRowInputs.length);
        var nnn = sampRows[sampRows.length - 1].name.split("_");
        newIndex = Number(nnn[nnn.length - 1]) + 1;
        newIndexStr = newIndex.toString();

        for (var i = 0; i < newRowInputs.length; i++) {
            n = newRowInputs[i].name.split("_")
            if (n[0] == "sampleID") {
                newRowInputs[i].value = "0";
            }
            newRowInputs[i].name = n[0] + "_" + n[1] + "_" + n[2] + "_" + n[3] + "_" + newIndexStr;
        }
        laneCell1.appendChild(newRow);
    }


    //Need to pass a whole data hierachy as an array

    //Organize diaplay as tables
    //One for each experiment contining two cells
    //	First cell is the experiment ID and info
    //	Second cell contains a table
    //		This table has two cells
    //		First cell is the Project details name...
    //		Second cell contains a table
    //			This table has two cells.
    //			First cell contains the lane data
    //			Second cell contains a table
    //				This table cont a row for each sample

    //	etc to organize the data
    // 	Come back later and add the ability to show and hide the data that is displayed

    var inForm = document.forms['upDateForm'];
    //var inForm=document.getElementById("upDateForm");
    var contentCell = document.getElementById("content");
    //contentCell.appendChild(inForm);
    //var numbLanes=8;
    for (var ex = 0; ex < exptArr.length; ex++)	//					Each experimenmt
    {
        var expTab = document.createElement("table");									//Make experiment table
        expTab.id = "expTab_" + ex.toString();
        expTab.border = 1;
        var expRow = document.createElement("tr");
        var expCell0 = document.createElement("td"); 				 //This is the cell that the experimental stuff goes in
        var expDataTab = document.createElement("table");
        for (var exData = 0; exData < exptArr[ex][0].length; exData++) {

            var rw = document.createElement("tr");
            var da = document.createElement("td");
            var cellText = document.createTextNode(exptArr[ex][0][exData][0]);
            //var cellText=document.createTextNode(exptArr[ex][0][exData][0]+"_"+ex.toString());
            da.appendChild(cellText);
            rw.appendChild(da);
            var da = document.createElement("td");
            if (exptArr[ex][0][exData][0] == "masterProjectID") {
                var cellIDtext = document.createTextNode(exptArr[ex][0][exData][1]);
                da.appendChild(cellIDtext);
                rw.appendChild(da);
                expDataTab.appendChild(rw);
                var cellIn = document.createElement('input');
                cellIn.type = 'hidden';
                cellIn.name = exptArr[ex][0][exData][0] + "_" + ex.toString();
                cellIn.value = exptArr[ex][0][exData][1];
                da.appendChild(cellIn);
                rw.appendChild(da);

            }
            else {
                var cellIn = document.createElement('input');
                cellIn.type = 'text';
                cellIn.name = exptArr[ex][0][exData][0] + "_" + ex.toString();
                cellIn.value = exptArr[ex][0][exData][1];
                da.appendChild(cellIn);
                rw.appendChild(da);
            }
            expDataTab.appendChild(rw);
        }

        expCell0.appendChild(expDataTab);
        expRow.appendChild(expCell0);
        //////////////////////////////////////////////////////////////////////////Run Cell1
        var expCell1 = document.createElement("td");  				//This is the cell that the rest go in

        for (var pr = 0; pr < exptArr[ex][1].length; pr++) {
            var projTab = document.createElement("table");
            projTab.id = "projTab_" + pr.toString();
            projTab.border = 1;
            var projRow = document.createElement("tr");
            ///////////////////////////////////////////////////////////////////////////////Pro Cell0
            var projCell0 = document.createElement("td"); 					//This is the cell that the project data goes in
            var projDataTab = document.createElement("table");
            for (var prData = 0; prData < exptArr[ex][1][pr][0].length; prData++) {
                var rw = document.createElement("tr");
                var da = document.createElement("td");
                var cellText = document.createTextNode(exptArr[ex][1][pr][0][prData][0]);
                //var cellText=document.createTextNode(exptArr[ex][1][pr][0][prData][0]+"_"+ex.toString()+"_"+pr.toString());
                da.appendChild(cellText);
                rw.appendChild(da);
                var da = document.createElement("td");
                if (exptArr[ex][1][pr][0][prData][0] == "seqProjectID") {
                    var cellIDtext = document.createTextNode(exptArr[ex][1][pr][0][prData][1]);
                    da.appendChild(cellIDtext);
                    rw.appendChild(da);
                    expDataTab.appendChild(rw);
                    var cellIn = document.createElement('input');
                    cellIn.type = 'hidden';
                    cellIn.name = exptArr[ex][1][pr][0][prData][0] + "_" + ex.toString() + "_" + pr.toString();
                    cellIn.value = exptArr[ex][1][pr][0][prData][1];
                    da.appendChild(cellIn);
                    rw.appendChild(da);
                }
                else {
                    var cellIn = document.createElement('input');
                    cellIn.type = 'text';
                    cellIn.name = exptArr[ex][1][pr][0][prData][0] + "_" + ex.toString() + "_" + pr.toString();
                    cellIn.value = exptArr[ex][1][pr][0][prData][1];
                    da.appendChild(cellIn);
                    rw.appendChild(da);
                }
                projDataTab.appendChild(rw);
            }
            projCell0.appendChild(projDataTab);
            projRow.appendChild(projCell0);
            //////////////////////////////////////////////////////////////////////////////////////Pro Cell1
            var projCell1 = document.createElement("td");  //This is the cell that the rest go in


            //var someText=document.createTextNode("Project name goes here"+pr.toString());

            //projCell.appendChild(someText);


            for (var la = 0; la < exptArr[ex][1][pr][1].length; la++) {
                var laneTab = document.createElement("table");
                laneTab.id = "laneTab+" + la.toString();
                laneTab.border = 1;
                var laneRow = document.createElement("tr");
                /////////////////////////////////////////////////////////////////////////////////////Lane Cell0
                var laneCell0 = document.createElement("td"); //This is the cell that the lane data goes in
                var laneDataTab = document.createElement("table");
                laneDataTab.border = 0;

                for (var laData = 0; laData < exptArr[ex][1][pr][1][la][0].length; laData++) {

                    var rw = document.createElement("tr");
                    var da = document.createElement("td");
                    var cellText = document.createTextNode(exptArr[ex][1][pr][1][la][0][laData][0]);
                    //var cellText=document.createTextNode(exptArr[ex][1][pr][1][la][0][laData][0]+"_"+ex.toString()+"_"+pr.toString()+"_"+la.toString());
                    da.appendChild(cellText);
                    rw.appendChild(da);
                    var da = document.createElement("td");
                    if (exptArr[ex][1][pr][1][la][0][laData][0] == "laneID") {
                        var cellIDtext = document.createTextNode(exptArr[ex][1][pr][1][la][0][laData][1]);
                        da.appendChild(cellIDtext);
                        rw.appendChild(da);
                        expDataTab.appendChild(rw);
                        var cellIn = document.createElement('input');
                        cellIn.type = 'hidden';
                        cellIn.name = exptArr[ex][1][pr][1][la][0][laData][0] + "_" + ex.toString() + "_" + pr.toString() + "_" + la.toString();
                        cellIn.value = exptArr[ex][1][pr][1][la][0][laData][1];
                        da.appendChild(cellIn);
                        rw.appendChild(da);
                    }
                    else {
                        var cellIn = document.createElement('input');
                        cellIn.type = 'text';
                        cellIn.name = exptArr[ex][1][pr][1][la][0][laData][0] + "_" + ex.toString() + "_" + pr.toString() + "_" + la.toString();
                        cellIn.value = exptArr[ex][1][pr][1][la][0][laData][1];
                        da.appendChild(cellIn);
                        rw.appendChild(da);
                    }
                    laneDataTab.appendChild(rw);
                }

                laneCell0.appendChild(laneDataTab);
                laneRow.appendChild(laneCell0);
                //////////////////////////////////////////////////////////////////////////////////////////////////Lane Cell1

                var laneCell1 = document.createElement("td");  //This contains just the samples. We put all the samples in lane cell1
                var lastSample;
                for (var sa = 0; sa < exptArr[ex][1][pr][1][la][1].length; sa++) {
                    var rw = document.createElement("tr");
                    rw.id = "sampleRow_" + ex.toString() + "_" + pr.toString() + "_" + la.toString() + "_" + sa.toString();
                    lastSample = ex.toString() + "_" + pr.toString() + "_" + la.toString() + "_" + sa.toString();

                    //alert(rw.id);
                    for (var saData = 0; saData < exptArr[ex][1][pr][1][la][1][sa].length; saData++) {
                        var da = document.createElement("td");
                        var cellText = document.createTextNode(exptArr[ex][1][pr][1][la][1][sa][saData][0]);
                        //var cellText=document.createTextNode(exptArr[ex][1][pr][1][la][1][sa][saData][0]+"_"+ex.toString()+"_"+pr.toString()+"_"+la.toString()+"_"+sa.toString());
                        da.appendChild(cellText);
                        rw.appendChild(da);
                        var da = document.createElement("td");
                        if (exptArr[ex][1][pr][1][la][1][sa][saData][0] == "sampleID") {
                            var cellIDtext = document.createTextNode(exptArr[ex][1][pr][1][la][1][sa][saData][1]);

                            da.appendChild(cellIDtext);
                            rw.appendChild(da);
                            //expDataTab.appendChild(rw);
                            var cellIn = document.createElement('input');
                            cellIn.id = "sampleID";
                            cellIn.type = 'hidden';
                            cellIn.name = exptArr[ex][1][pr][1][la][1][sa][saData][0] + "_" + ex.toString() + "_" + pr.toString() + "_" + la.toString() + "_" + sa.toString();
                            cellIn.value = exptArr[ex][1][pr][1][la][1][sa][saData][1];
                            da.appendChild(cellIn);
                            rw.appendChild(da);
                        }
                        else {
                            var cellIn = document.createElement('input');
                            cellIn.type = 'text';
                            cellIn.name = exptArr[ex][1][pr][1][la][1][sa][saData][0] + "_" + ex.toString() + "_" + pr.toString() + "_" + la.toString() + "_" + sa.toString();
                            cellIn.value = exptArr[ex][1][pr][1][la][1][sa][saData][1];
                            da.appendChild(cellIn);
                            //da.appendChild(cellText);
                            rw.appendChild(da);
                        }
                    }

                    laneCell1.appendChild(rw);
                }
                laneCell1.id = "laneCell_" + lastSample;
                var rw = document.createElement("tr");
                var da = document.createElement("td");
                var addSampButton = document.createElement('input');
                addSampButton.type = 'button';
                addSampButton.value = 'Add sample';
                addSampButton.name = "sampleButton_" + lastSample;
                addSampButton.id = lastSample;
                //alert(addSampButton.name);
                //use the "this" keywork to pass attributes of the current element.
                //
                addSampButton.addEventListener('click', function () {
                    addSample(this.id);
                }, false);

                da.appendChild(addSampButton);
                rw.appendChild(da);
                laneCell1.appendChild(rw);


                laneRow.appendChild(laneCell1);
                laneTab.appendChild(laneRow);

                projCell1.appendChild(laneTab);
            }

            projRow.appendChild(projCell1);
            projTab.appendChild(projRow);

            expCell1.appendChild(projTab);
            //var someText.createTextNode(pr.toString());
            //projCell.appendChild(someText);


        }


        expRow.appendChild(expCell1);
        expTab.appendChild(expRow);

        inForm.appendChild(expTab);
        //contentCell.appendChild(inForm);
        //contentCell.appendChild(expTab);
    }

}
