/**
 * Created by jon on 24/11/15.
 */

function addSelectionGroupBehaviours() {
    $('.selection').change(function() {
        var selected = $(this).find(':selected').text();
        $(".description").hide();
        $("[id='" + selected + "'").show();
        selected == "None" ? $('#toggle').show() : $('#toggle').hide();
        window.dispatchEvent(new Event('resize'));
        $(window).trigger('resize');
    });
}

$(document).ready(function() {
    addSelectionGroupBehaviours();
    $('.selection').val("None")
});