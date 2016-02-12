/**
 * Created by jon on 10/02/16.
 */

function addTooltipBehaviours() {
    $('input').tooltip({
        placement: "right",
        trigger: "hover"
    });
}

$(document).ready(function() {
    addTooltipBehaviours();
});