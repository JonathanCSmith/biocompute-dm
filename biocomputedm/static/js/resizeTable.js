/**
 * Created by jon on 15/02/16.
 */

function addResizeBehaviours() {
    $(window).on("resize", resize());
}

function resize() {
    var size =$(".resize").width();
    if (size > 200) {

        $(".resize").height(200);
        $(".resize").width(200);
    }

    else {
        $(".resize").height(size);
    }
}

$(document).ready(function() {
    resize();
    addResizeBehaviours();
});