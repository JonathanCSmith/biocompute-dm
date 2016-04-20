/**
 * Created by jon on 15/02/16.
 */

function addResizeBehaviours() {
    $(window).on("resize", resizer);
}

function resizer() {
    var size =$(".resize:visible").width();
    if (size > 200) {
        $(".resize:visible").height(200);
        $(".resize:visible").width(200);
        $(".resize_header:visible").width(200);
    }

    else {
        $(".resize:visible").height(size);
        $(".resize_header:visible").width(size);
    }

    return;
}

$(document).ready(function() {
    resizer();
    addResizeBehaviours();
});