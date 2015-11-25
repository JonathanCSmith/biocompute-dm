/**
 * Created by jon on 24/11/15.
 */

function addSelectionGroupBehaviours() {
    $('.selection').change(function(){
        var selected = $(this).find(':selected').text();
        //alert(selected);
        $(".description").hide();
        $('#' + selected).show();
    });
}

$(document).ready(function() {
   addSelectionGroupBehaviours();
});