

function experimentTable(exptArr)
	{

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
        function addSample(laneNo)
                {
                var frm=document.getElementById("newSeqRun");
                frm.action="/cgi-bin/coreInSys/updateAddSample";
                var laneAddTo=document.createElement("input");
                laneAddTo.type='hidden';
                laneAddTo.name="laneAddTo";
                laneAddTo.value=laneNo;
                frm.appendChild(laneAddTo);

		var nTA=document.getElementById("numSamps_"+laneNo);
		numbToAdd=document.createElement("input");
		numbToAdd.type='hidden';
		numbToAdd.name="numbToAdd";
		numbToAdd.value=nTA.value;
		frm.appendChild(numbToAdd);
                frm.submit();

                }
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	function remSample(sampNo)
                {
                var frm=document.getElementById("newSeqRun");
                frm.action="/cgi-bin/coreInSys/updateRemSample";
                var sampToRem=document.createElement("input");
                sampToRem.type='hidden';
                sampToRem.name="sampToRem";
                sampToRem.value=sampNo;
                frm.appendChild(sampToRem);

                frm.submit();

                }




	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	function laneCell1(exptArr,ex,pr,la)
		{
		var laneCell1=document.createElement("td");  //This contains just the samples. We put all the samples in lane cell1
		var lastSample;
                for (var sa=0;sa<exptArr[ex][1][pr][1][la][1].length;sa++)
               		{ 
                        lastSample=sa.toString();
                        var rw=document.createElement("tr");
			lastSample=ex.toString()+"_"+pr.toString()+"_"+la.toString()+"_"+sa.toString();
			rw.id="sampleRow_"+lastSample;

                        for (var saData=0;saData<exptArr[ex][1][pr][1][la][1][sa].length;saData++)
                        	{

                                var da=document.createElement("td");
                                var cellText=document.createTextNode(exptArr[ex][1][pr][1][la][1][sa][saData][0]);

				da.appendChild(cellText);
                                rw.appendChild(da);
                                var da=document.createElement("td");
                                if (exptArr[ex][1][pr][1][la][1][sa][saData][0]=="sampleID")
                                	{
                                        var cellIDtext=document.createTextNode(exptArr[ex][1][pr][1][la][1][sa][saData][1]);
                                        da.appendChild(cellIDtext);
                                        rw.appendChild(da);
					var cellIn=document.createElement('input');
                                        cellIn.type='hidden';
                                        cellIn.name=exptArr[ex][1][pr][1][la][1][sa][saData][0]+"_"+ex.toString()+"_"+pr.toString()+"_"+la.toString()+"_"+sa.toString();
                                        cellIn.value=exptArr[ex][1][pr][1][la][1][sa][saData][1];
                                        da.appendChild(cellIn);
                                        rw.appendChild(da);
                                        }
				else
                                        {
                                        var cellIn=document.createElement('input');
                                        cellIn.type='text';
                                        cellIn.name=exptArr[ex][1][pr][1][la][1][sa][saData][0]+"_"+ex.toString()+"_"+pr.toString()+"_"+la.toString()+"_"+sa.toString();
                                        cellIn.value=exptArr[ex][1][pr][1][la][1][sa][saData][1];
                                        da.appendChild(cellIn);
					rw.appendChild(da);
                                        }
				}
			//add a remove sample button if there are multiple samples in this lane
			if (exptArr[ex][1][pr][1][la][1].length>1)
                        	{
                        	//alert(exptArr[ex][1][pr][1][la][1].length);
				var da=document.createElement("td");
                		var remSampButton=document.createElement('input');
                		remSampButton.type='button';
                		remSampButton.value='Remove sample';
                		remSampButton.name="remSampleButton_"+lastSample;
                                remSampButton.id=lastSample;
                                remSampButton.addEventListener('click',function() {remSample(this.id);}, false);
           
                                da.appendChild(remSampButton);
                                rw.appendChild(da);
	
                        	}



			laneCell1.appendChild(rw);

                	}

		laneCell1.id="laneCell_"+lastSample;
                var rw=document.createElement("tr");
                var da=document.createElement("td");
                var selectNumSamps=document.createElement('select');

                //addSampButton.type='select';
		for (var samNo=1;samNo<17;samNo++)
			{	
			var option=document.createElement("option");
			numbOpt=samNo.toString();
			option.text=numbOpt;	
			option.value=numbOpt;
			selectNumSamps.add(option,null);
			}
		selectNumSamps.id="numSamps_"+lastSample;
		//alert(selectNumSamps.id);
		da.appendChild(selectNumSamps);
		var addSampButton=document.createElement('input');
		addSampButton.type='button';
                addSampButton.value='Add sample';
                addSampButton.name="sampleButton_"+lastSample;
		//addSampButton.name="sampleButton";
                addSampButton.id=lastSample;

            	addSampButton.addEventListener('click',function() {addSample(this.id);}, false);

                da.appendChild(addSampButton);
                rw.appendChild(da);
                laneCell1.appendChild(rw);

		return(laneCell1);

		}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	function makeLaneCell0(exptArr,ex,pr,la,lastLane)
		{
		var laneCell0=document.createElement("td"); //This is the cell that the lane data goes in
                var laneDataTab=document.createElement("table");

                laneDataTab.border=0;
		var laneID=0;
		for (var laData=0;laData<exptArr[ex][1][pr][1][la][0].length;laData++)
			{
                        var rw=document.createElement("tr");
                        var da=document.createElement("td");

                        var cellText=document.createTextNode(exptArr[ex][1][pr][1][la][0][laData][0]);
                        da.appendChild(cellText);
                        rw.appendChild(da);
                        var da=document.createElement("td");
                        if (exptArr[ex][1][pr][1][la][0][laData][0]=="laneID")
                        	{
				laneID=exptArr[ex][1][pr][1][la][0][laData][1];
                                var cellIDtext=document.createTextNode(exptArr[ex][1][pr][1][la][0][laData][1]);
                                da.appendChild(cellIDtext);
                                rw.appendChild(da);
                                expDataTab.appendChild(rw);
                                var cellIn=document.createElement('input');
                                cellIn.type='hidden';
                                cellIn.name=exptArr[ex][1][pr][1][la][0][laData][0]+"_"+ex.toString()+"_"+pr.toString()+"_"+la.toString();
                                cellIn.value=exptArr[ex][1][pr][1][la][0][laData][1];
                                da.appendChild(cellIn);
                                rw.appendChild(da);
				}
			else
                        	{
                                var cellIn=document.createElement('input');
                                cellIn.type='text';
                                cellIn.name=exptArr[ex][1][pr][1][la][0][laData][0]+"_"+ex.toString()+"_"+pr.toString()+"_"+la.toString();
                                cellIn.value=exptArr[ex][1][pr][1][la][0][laData][1];
                                da.appendChild(cellIn);
                                rw.appendChild(da);
                                }
			laneDataTab.appendChild(rw);
                        }
		if (exptArr[ex][1][pr][1].length>1)
                        {
			var rw1=document.createElement("tr");
                        var da1=document.createElement("td");
                        var remLaneButton=document.createElement('input');
                        remLaneButton.type='button';
                        remLaneButton.value='Remove lane';
                        remLaneButton.name="remLaneButton_"+lastLane;
                        remLaneButton.id=lastLane;
                        //remLaneButton.addEventListener('click',function() {remLane(this.id,laneID);}, false);
			remLaneButton.addEventListener('click',function() {remLane(this.id,laneID);}, false);
                        da1.appendChild(remLaneButton);
                        rw1.appendChild(da1);
			laneDataTab.appendChild(rw1);
                        }




		laneCell0.appendChild(laneDataTab);

		return(laneCell0);

		}




	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	function addLane(projNo)
		{
		//get the form
		var frm=document.getElementById("newSeqRun");
		//change the form's action
		frm.action="/cgi-bin/coreInSys/updateAddLane";

		//save the project number that we want to append the lane to
		var projAddTo=document.createElement("input");
		projAddTo.type='hidden';
		projAddTo.name="projAddTo";
		projAddTo.value=projNo;
		//alert(projAddTo.value);
		frm.appendChild(projAddTo);
		//submit the form
		frm.submit();

		}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//function remLane(laneNo,laneID)
	function remLane(laneNo)
                {
		//alert(laneNo);
		//alert(laneID);
                var frm=document.getElementById("newSeqRun");
                frm.action="/cgi-bin/coreInSys/updateRemLane";
                var laneToRem=document.createElement("input");
                laneToRem.type='hidden';
                laneToRem.name="laneToRem";
                laneToRem.value=laneNo;
                frm.appendChild(laneToRem);
                //frm.submit();
		//var rmlLen=remLanes.length;
		//remLanes[rmlLen+1]=laneID;
		frm.submit();
                }
	


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	function addSeqProj()
		{
		var frm=document.getElementById("newSeqRun");
		frm.action="/cgi-bin/coreInSys/updateAddSeqPro";		

		frm.submit();

		}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        function remSeqProj(projNo)
                {
                var frm=document.getElementById("newSeqRun");
                frm.action="/cgi-bin/coreInSys/updateRemSeqProj";
                var projToRem=document.createElement("input");
                projToRem.type='hidden';
                projToRem.name="seqProjToRem";
                projToRem.value=projNo;
                frm.appendChild(projToRem);
                frm.submit();
                }





	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	function projCell1(exptArr,ex,pr)
		{
		var projCell1=document.createElement("td");  //This is the cell that the rest go in
                var lastLane;
		var laneCell0;		
		var rw=document.createElement("tr");
                var da=document.createElement("td");
                for (var la=0;la<exptArr[ex][1][pr][1].length;la++)
                	{
                        var laneTab=document.createElement("table");
			lastLane=ex.toString()+"_"+pr.toString()+"_"+la.toString();
                        //projCell1.id="laneDataTab_"+lastLane;
                        laneTab.border=1;
                        var laneRow=document.createElement("tr");

			laneCell0=makeLaneCell0(exptArr,ex,pr,la,lastLane);

                        laneRow.appendChild(laneCell0);
 			var LC1=laneCell1(exptArr,ex,pr,la);
                        laneRow.appendChild(LC1);
                        laneTab.appendChild(laneRow);
                        //projCell1.appendChild(laneTab);
			da.appendChild(laneTab);
			rw.appendChild(da);
			rw.id="laneDataTab_"+lastLane;

			projCell1.appendChild(rw);	


	
                        }
		projCell1.id="projCell1_"+lastLane;
                var rw=document.createElement("tr");
                var da=document.createElement("td");
                var addLaneButton=document.createElement('input');
                addLaneButton.type='button';
                addLaneButton.value='Add lane';
		addLaneButton.name="laneButton_"+lastLane;

		//var prBit=lastLane.split("_");
		//addLaneButton.id=prBit[1];
 		addLaneButton.id=lastLane;
		//get the form

                addLaneButton.addEventListener('click',function() {addLane(this.id);}, false);

                da.appendChild(addLaneButton);
                rw.appendChild(da);
                laneCell0.appendChild(rw);

		return(projCell1);
		}




	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	function makeProjectCell0(exptArr,ex,pr,lastProj)
		{
		var projCell0=document.createElement("td");
                var projDataTab=document.createElement("table");
                for (var prData=0;prData<exptArr[ex][1][pr][0].length;prData++)
                	{
                        var rw=document.createElement("tr");
                        var da=document.createElement("td");
                        var cellText=document.createTextNode(exptArr[ex][1][pr][0][prData][0]);
                        da.appendChild(cellText);
                        rw.appendChild(da);
                        var da=document.createElement("td");
                        if (exptArr[ex][1][pr][0][prData][0]=="seqProjectID")
                        	{
                                var cellIDtext=document.createTextNode(exptArr[ex][1][pr][0][prData][1]);
                                da.appendChild(cellIDtext);
                                rw.appendChild(da);
                                expDataTab.appendChild(rw);
                                var cellIn=document.createElement('input');
                                cellIn.type='hidden';
                                cellIn.name=exptArr[ex][1][pr][0][prData][0]+"_"+ex.toString()+"_"+pr.toString();
                                cellIn.value=exptArr[ex][1][pr][0][prData][1];
                                da.appendChild(cellIn);
                                rw.appendChild(da);
                                }
			else
                        	{
                                var cellIn=document.createElement('input');
                                cellIn.type='text';
                                cellIn.name=exptArr[ex][1][pr][0][prData][0]+"_"+ex.toString()+"_"+pr.toString();
                                cellIn.value=exptArr[ex][1][pr][0][prData][1];
                                da.appendChild(cellIn);
                                rw.appendChild(da);
                                }
			projDataTab.appendChild(rw);
                        }


		 if (exptArr[ex][1].length>1)
                        {
                        var rw1=document.createElement("tr");
                        var da1=document.createElement("td");
                        var remProjButton=document.createElement('input');
                        remProjButton.type='button';
                        remProjButton.value='Remove';
                        remProjButton.name="remLaneButton_"+lastProj;
                        remProjButton.id=lastProj;
                        remProjButton.addEventListener('click',function() {remSeqProj(this.id);}, false);

                        da1.appendChild(remProjButton);
                        rw1.appendChild(da1);
                        projDataTab.appendChild(rw1);
                        }




		projCell0.appendChild(projDataTab);
		return (projCell0);	
		
		}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	function expCell1(exptArr,ex)
		{
		var expCell1=document.createElement("td");
		var lastProj;
                for (var pr=0;pr<exptArr[ex][1].length;pr++)
                        {
                        var projTab=document.createElement("table");
			lastLane=ex.toString()+"_"+pr.toString();
                        projTab.id="projTab_"+pr.toString();
                        projTab.border=1;
                        var projRow=document.createElement("tr");
			var projCell0=makeProjectCell0(exptArr,ex,pr,lastLane);

                        projRow.appendChild(projCell0);
                        var PC1=projCell1(exptArr,ex,pr);
                        projRow.appendChild(PC1);
                        projTab.appendChild(projRow);

                        expCell1.appendChild(projTab);

			}
		expCell1.id="expCell1_"+lastLane;
                var rw=document.createElement("tr");
                var da=document.createElement("td");
                var addSeqProButton=document.createElement('input');
                addSeqProButton.type='button';
                addSeqProButton.value='Add seqProj';
                addSeqProButton.name="seqProButton_"+lastLane;

                var prBit=lastLane.split("_");
                addSeqProButton.id=prBit[1];
                addSeqProButton.id=lastLane;
                //get the form
                
                addSeqProButton.addEventListener('click',function() {addSeqProj(this.id);}, false);
           
                da.appendChild(addSeqProButton);
                rw.appendChild(da);
                projCell0.appendChild(rw);





		return(expCell1);
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

	var inForm=document.forms['newSeqRun'];
/*
for (var prData=0;prData<exptArr[ex][1][pr][0].length;prData++)
                        {
                        var rw=document.createElement("tr");
                        var da=document.createElement("td");
                        var cellText=document.createTextNode(exptArr[ex][1][pr][0][prData][0]);
                        da.appendChild(cellText);
                        rw.appendChild(da);
                        var da=document.createElement("td");
                        if (exptArr[ex][1][pr][0][prData][0]=="seqProjectID")
                                {
                                var cellIDtext=document.createTextNode(exptArr[ex][1][pr][0][prData][1]);
                                da.appendChild(cellIDtext);
                                rw.appendChild(da);
                                expDataTab.appendChild(rw);
                                var cellIn=document.createElement('input');
                                cellIn.type='hidden';
                                cellIn.name=exptArr[ex][1][pr][0][prData][0]+"_"+ex.toString()+"_"+pr.toString();
                                cellIn.value=exptArr[ex][1][pr][0][prData][1];
                                da.appendChild(cellIn);
                                rw.appendChild(da);
                                }
                        else
                                {
                                var cellIn=document.createElement('input');
                                cellIn.type='text';
                                cellIn.name=exptArr[ex][1][pr][0][prData][0]+"_"+ex.toString()+"_"+pr.toString();
                                cellIn.value=exptArr[ex][1][pr][0][prData][1];
                                da.appendChild(cellIn);
                                rw.appendChild(da);
                                }
                        projDataTab.appendChild(rw);
                        }
*/
	



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
			//alert (exptArr[ex][0][exData][0]);
			if (exptArr[ex][0][exData][0]=="seqRunID")

				{
				var cellIDtext=document.createTextNode(exptArr[ex][0][exData][1]);
                                da.appendChild(cellIDtext);
                                rw.appendChild(da);
                                expDataTab.appendChild(rw);
                                var cellIn=document.createElement('input');
                                cellIn.type='hidden';
                                cellIn.name=exptArr[ex][0][exData][0]+"_"+ex.toString();
                                cellIn.value=exptArr[ex][0][exData][1];
                                da.appendChild(cellIn);
                                rw.appendChild(da);
                                }
                        else
                                {
				var cellIn=document.createElement('input');
				cellIn.type='text';
				cellIn.name=exptArr[ex][0][exData][0]+"_"+ex.toString();
				cellIn.value=exptArr[ex][0][exData][1];
				//inForm.appendChild(cellIn);
                        	da.appendChild(cellIn);
                        	rw.appendChild(da);
				}
			expDataTab.appendChild(rw);
				
			}

		expCell0.appendChild(expDataTab);
		expRow.appendChild(expCell0);
		//////////////////////////////////////////////////////////////////////////Run Cell1

		var EC1=expCell1(exptArr,ex);   		
		

		expRow.appendChild(EC1);
		expTab.appendChild(expRow);
	
		inForm.appendChild(expTab);


		
		}	


	
	}
