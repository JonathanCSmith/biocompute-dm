/**
 * Created by jon on 24/11/15.
 */

function addSelectionGroupBehaviours() {
    $('.selection').change(function(){
        var selected = $(this).find(':selected').text();
        $(".description").hide();
        $("[id='" + selected + "'").show();
    });
}

$(document).ready(function() {
   addSelectionGroupBehaviours();
});