/**
 * Created by jon on 10/02/16.
 */

function addTooltipBehaviours() {
    $(document).tooltip({
        placement: "right",
        trigger: "hover"
    });
}

$(document).ready(function() {
    addTooltipBehaviours();
});