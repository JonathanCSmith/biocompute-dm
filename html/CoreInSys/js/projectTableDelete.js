


function experimentTable(exptArr)
	{
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

	
	var contentCell=document.getElementById("content");
	//var numbLanes=8;
	for (var ex=0;ex<exptArr.length;ex++)	//					Each experimenmt
		{
		var expTab=document.createElement("table");									//Make experiment table
		expTab.id="expTab_"+ex.toString();
		expTab.border=1;
		var expRow=document.createElement("tr");
		///////////////////////////////////////////////////////////////////////////Run Cell0
		var expCell0=document.createElement("td"); 				 //This is the cell that the experimental stuff goes in
		var expDataTab=document.createElement("table");
		for (var exData=0;exData<exptArr[ex][0].length;exData++)
			{

			if (exptArr[ex][0][exData][0]=="masterProjectID")
                        	{
                                var rw=document.createElement("tr");
                                var delcell=document.createElement("td");
                                var a = document.createElement('a');
                                var linkText = document.createTextNode("Remove master project from DB");
                                a.appendChild(linkText);
                                a.title = "delete masterProject";
                                a.href = "delMasterProj?masterProjectID="+exptArr[ex][0][exData][1];
                                delcell.appendChild(a);
                                rw.appendChild(delcell);
                                expDataTab.appendChild(rw);
                                }

			var rw=document.createElement("tr");
			var da=document.createElement("td");
			var cellText=document.createTextNode(exptArr[ex][0][exData][0]);
			da.appendChild(cellText);
			rw.appendChild(da);
			var da=document.createElement("td");
                        var cellText=document.createTextNode(exptArr[ex][0][exData][1]);
                        da.appendChild(cellText);
                        rw.appendChild(da);
			expDataTab.appendChild(rw);
			}

		expCell0.appendChild(expDataTab);
		expRow.appendChild(expCell0);
		//////////////////////////////////////////////////////////////////////////Run Cell1
		var expCell1=document.createElement("td");  				//This is the cell that the rest go in
		
		for (var pr=0;pr<exptArr[ex][1].length;pr++)
			{
			var projTab=document.createElement("table");
			projTab.id="projTab_"+pr.toString();
			projTab.border=1;
			var projRow=document.createElement("tr");
			///////////////////////////////////////////////////////////////////////////////Pro Cell0
			var projCell0=document.createElement("td"); 					//This is the cell that the project data goes in
			var projDataTab=document.createElement("table");
			for (var prData=0;prData<exptArr[ex][1][pr][0].length;prData++)
                        	{
				if (exptArr[ex][1][pr][0][prData][0]=="seqProjectID")
					{
					var rw=document.createElement("tr");
                                        var delcell=document.createElement("td");
					var a = document.createElement('a');
                                        var linkText = document.createTextNode("Remove seq project from DB");
                                        a.appendChild(linkText);
                                        a.title = "delete seqProject";
                                        a.href = "removeSeqProj?seqProjectID="+exptArr[ex][1][pr][0][prData][1];
                                        delcell.appendChild(a);
					rw.appendChild(delcell);
                                        projDataTab.appendChild(rw);
					}

                        	var rw=document.createElement("tr");
                        	var da=document.createElement("td");
                        	var cellText=document.createTextNode(exptArr[ex][1][pr][0][prData][0]);
                        	da.appendChild(cellText);
                        	rw.appendChild(da);
                        	var da=document.createElement("td");
                        	var cellText=document.createTextNode(exptArr[ex][1][pr][0][prData][1]);
                        	da.appendChild(cellText);
                        	rw.appendChild(da);
                        	projDataTab.appendChild(rw);
                        	}
			projCell0.appendChild(projDataTab);
			projRow.appendChild(projCell0);
			//////////////////////////////////////////////////////////////////////////////////////Pro Cell1
			var projCell1=document.createElement("td");  //This is the cell that the rest go in



			//var someText=document.createTextNode("Project name goes here"+pr.toString());

                        //projCell.appendChild(someText);
			

			for (var la=0;la<exptArr[ex][1][pr][1].length;la++)
				{
				var laneTab=document.createElement("table");
				laneTab.id="laneTab+"+la.toString();
				laneTab.border=1;
				var laneRow=document.createElement("tr");
				/////////////////////////////////////////////////////////////////////////////////////Lane Cell0
				var laneCell0=document.createElement("td"); //This is the cell that the lane data goes in
				var laneDataTab=document.createElement("table");
				laneDataTab.border=0;
				
				for (var laData=0;laData<exptArr[ex][1][pr][1][la][0].length;laData++)
					{

					if (exptArr[ex][1][pr][1][la][0][laData][0]=="laneID")
						{	
						var rw=document.createElement("tr");
						var delcell=document.createElement("td");
						//var delcellTest=document.createTextNode("DELETE");
						var a = document.createElement('a');
						var linkText = document.createTextNode("Remove lane from DB");
						a.appendChild(linkText);
						a.title = "delete lane";
						a.href = "removeLane?laneID="+exptArr[ex][1][pr][1][la][0][laData][1];
						delcell.appendChild(a);

						//delcell.appendChild(delcellTest);
						rw.appendChild(delcell);
						laneDataTab.appendChild(rw);
						}
						
					var rw=document.createElement("tr");
                                	var da=document.createElement("td");
                                	var cellText=document.createTextNode(exptArr[ex][1][pr][1][la][0][laData][0]);
				
                                	da.appendChild(cellText);
                                	rw.appendChild(da);
                                	var da=document.createElement("td");
                                	var cellText=document.createTextNode(exptArr[ex][1][pr][1][la][0][laData][1]);
                                	da.appendChild(cellText);
                                	rw.appendChild(da);
                                	laneDataTab.appendChild(rw);
					}

				laneCell0.appendChild(laneDataTab);
				laneRow.appendChild(laneCell0);
				//////////////////////////////////////////////////////////////////////////////////////////////////Lane Cell1	
				
				var laneCell1=document.createElement("td");  //This contains just the samples. We put all the samples in lane cell1


				for (var sa=0;sa<exptArr[ex][1][pr][1][la][1].length;sa++)
					{
					var rw=document.createElement("tr");
					for (var saData=0;saData<exptArr[ex][1][pr][1][la][1][sa].length;saData++)
						{
						if (exptArr[ex][1][pr][1][la][1][sa][saData][0]=="sampleID")
							{
							var delcell=document.createElement("td");
							var a = document.createElement('a');
							var linkText = document.createTextNode("Remove sample from DB");
							a.appendChild(linkText);
                                                	a.title = "delete sample";
							a.href = "removeSample?sampleID="+exptArr[ex][1][pr][1][la][1][sa][saData][1];
							delcell.appendChild(a);
							rw.appendChild(delcell);
							}

                                        	var da=document.createElement("td");
                                        	var cellText=document.createTextNode(exptArr[ex][1][pr][1][la][1][sa][saData][0]);
                                        	da.appendChild(cellText);
                                        	rw.appendChild(da);
                                        	var da=document.createElement("td");
                                        	var cellText=document.createTextNode(exptArr[ex][1][pr][1][la][1][sa][saData][1]);
                                        	da.appendChild(cellText);
                                        	rw.appendChild(da);
						}
					laneCell1.appendChild(rw);
					}
				
				
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
		
		contentCell.appendChild(expTab);
		}	
	
	}



