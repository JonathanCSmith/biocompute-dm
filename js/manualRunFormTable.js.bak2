

function experimentTable(exptArr)
	{



	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
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

	var inForm=document.forms['newSeqRun'];

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


			var rw=document.createElement("tr");
			var da=document.createElement("td");
			var cellText=document.createTextNode(exptArr[ex][0][exData][0]);
			//var cellText=document.createTextNode(exptArr[ex][0][exData][0]+"_"+ex.toString());
			da.appendChild(cellText);
			rw.appendChild(da);
			var da=document.createElement("td");
			if (exptArr[ex][0][exData][0]=="seqRunID")
				{
				var cellIn=document.createElement('p');
				}
			else
				{
				var cellIn=document.createElement('input');
				}
			cellIn.type='text';
			cellIn.name=exptArr[ex][0][exData][0]+"_"+ex.toString();
			cellIn.value=exptArr[ex][0][exData][1];
			inForm.appendChild(cellIn);
                        da.appendChild(cellIn);
                        rw.appendChild(da);
			expDataTab.appendChild(rw);
			}

		expCell0.appendChild(expDataTab);
		expRow.appendChild(expCell0);
		//////////////////////////////////////////////////////////////////////////Run Cell1


		expTab.appendChild(expRow);
	
		inForm.appendChild(expTab);


		
		}	
	
	}
