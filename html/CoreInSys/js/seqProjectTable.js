


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

	
	//var contentCell=document.getElementById("content");
	var runDetailsShown=document.getElementById("runDetailsShown");
		
	//var numbLanes=8;
		
		for (var pr=0;pr<exptArr.length;pr++)
			{
			var projTab=document.createElement("table");
			projTab.id="projTab_"+pr.toString();
			projTab.border=1;
			var projRow=document.createElement("tr");
			///////////////////////////////////////////////////////////////////////////////Pro Cell0
			var projCell0=document.createElement("td"); 					//This is the cell that the project data goes in
			projCell0.style.verticalAlign="top";
			var projDataTab=document.createElement("table");
			for (var prData=0;prData<exptArr[pr][0].length;prData++)
                        	{
                        	var rw=document.createElement("tr");
                        	var da=document.createElement("td");
                        	var cellText=document.createTextNode(exptArr[pr][0][prData][0]);
                        	da.appendChild(cellText);
                        	rw.appendChild(da);
                        	var da=document.createElement("td");
                        	var cellText=document.createTextNode(exptArr[pr][0][prData][1]);
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
			

			for (var la=0;la<exptArr[pr][1].length;la++)
				{
				var laneTab=document.createElement("table");
				laneTab.id="laneTab+"+la.toString();
				laneTab.border=1;
				var laneRow=document.createElement("tr");
				/////////////////////////////////////////////////////////////////////////////////////Lane Cell0
				var laneCell0=document.createElement("td"); //This is the cell that the lane data goes in
				laneCell0.style.verticalAlign="top";
				var laneDataTab=document.createElement("table");
				laneDataTab.border=0;
				
				for (var laData=0;laData<exptArr[pr][1][la][0].length;laData++)
					{

					var rw=document.createElement("tr");
                                	var da=document.createElement("td");
                                	var cellText=document.createTextNode(exptArr[pr][1][la][0][laData][0]);
				
                                	da.appendChild(cellText);
                                	rw.appendChild(da);
                                	var da=document.createElement("td");
                                	var cellText=document.createTextNode(exptArr[pr][1][la][0][laData][1]);
                                	da.appendChild(cellText);
                                	rw.appendChild(da);
                                	laneDataTab.appendChild(rw);
					}

				laneCell0.appendChild(laneDataTab);
				laneRow.appendChild(laneCell0);
				//////////////////////////////////////////////////////////////////////////////////////////////////Lane Cell1	
				
				var laneCell1=document.createElement("td");  //This contains just the samples. We put all the samples in lane cell1


				for (var sa=0;sa<exptArr[pr][1][la][1].length;sa++)
					{
					var rw=document.createElement("tr");
					for (var saData=0;saData<exptArr[pr][1][la][1][sa].length;saData++)
						{

                                        	var da=document.createElement("td");
                                        	var cellText=document.createTextNode(exptArr[pr][1][la][1][sa][saData][0]);
                                        	da.appendChild(cellText);
                                        	rw.appendChild(da);
                                        	var da=document.createElement("td");
                                        	var cellText=document.createTextNode(exptArr[pr][1][la][1][sa][saData][1]);
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



		
		runDetailsShown.appendChild(projTab);
		}	
	
	}
