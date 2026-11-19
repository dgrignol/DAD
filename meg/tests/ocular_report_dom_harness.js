/*
 * Hermetic source-level UI contract test; no browser, URL or network is used.
 *
 * Usage: node ocular_report_dom_harness.js VIEWER_JS SYNTHETIC_PAYLOAD_JSON
 * The Python test supplies small synthetic traces. A minimal DOM/canvas stub
 * records drawing operations and invokes the viewer's own event handlers.
 * This checks navigation state, decoding and controls, not browser rendering.
 */
"use strict";
const fs = require("fs"), vm = require("vm"), assert = require("assert");
const payload = JSON.parse(fs.readFileSync(process.argv[3], "utf8"));
const draws = [], handlers = {}, elements = {};
const context = new Proxy({}, {get: function (target, key) {
  if (!(key in target)) target[key] = function () { draws.push([key].concat(Array.from(arguments))); };
  return target[key];
}});
function element(id) {
  if (!elements[id]) elements[id] = {id: id, value: "", textContent: "", children: [],
    hidden: false, disabled: false, style: {}, parentElement: {clientWidth: 1000},
    clientWidth: 1000, tagName: "DIV",
    appendChild: function (child) { this.children.push(child); if (!this.value && child.value) this.value = child.value; },
    replaceChildren: function () { this.children = []; this.value = ""; },
    getContext: function () { return context; }, setPointerCapture: function () {},
    addEventListener: function (kind, fn) { this[kind] = fn; }};
  return elements[id];
}
element("report-data").textContent = JSON.stringify(payload);
const document = {getElementById: element, createElement: function (tag) {
  return {tagName: tag.toUpperCase(), children: [], appendChild: function (child) { this.children.push(child); }};
}, documentElement: {dataset: {}}, activeElement: {tagName: "BODY"}};
const window = {devicePixelRatio: 1, addEventListener: function (kind, fn) { handlers[kind] = fn; }};
vm.runInNewContext(fs.readFileSync(process.argv[2], "utf8"), {
  document: document, window: window, Map: Map, Float32Array: Float32Array, Uint16Array: Uint16Array,
  Uint8Array: Uint8Array, atob: function (text) { return Buffer.from(text, "base64").toString("binary"); }
});
assert.strictEqual(document.documentElement.dataset.reportReady, "true");
const ready = payload.runs.filter(function (r) { return r.status === "ready"; });
if (!ready.length) {
  assert.strictEqual(element("review").hidden, true);
  assert.strictEqual(element("empty").hidden, false);
} else {
  assert.strictEqual(element("source").textContent, ready[0].source_filename);
  assert.ok(draws.some(function (r) { return r[0] === "lineTo"; }));
  element("window").value = "2"; element("window").onchange();
  element("jump").value = "2.5"; element("jump").onchange();
  assert.ok(element("time-label").textContent.startsWith("2.50 to 4.50"));
  element("forward").onclick();
  assert.ok(Number(element("time-slider").value) > 2.5);
  element("back").onclick();
  element("zoom-in").onclick();
  assert.ok(element("time-label").textContent.includes("3.00 to 4.00"));
  if (ready.length > 1) {
    element("next-run").onclick();
    assert.strictEqual(element("source").textContent, ready[1].source_filename);
    assert.strictEqual(element("previous-run").disabled, false);
    assert.strictEqual(element("next-run").disabled, true);
    element("previous-run").onclick();
    assert.strictEqual(element("source").textContent, ready[0].source_filename);
  }
  element("window").value = "2"; element("window").onchange();
  const before = Number(element("time-slider").value);
  element("traces").onpointerdown({clientX: 400, pointerId: 1});
  element("traces").onpointermove({clientX: 200});
  element("traces").onpointerup();
  assert.ok(Number(element("time-slider").value) > before);
  handlers.keydown({key: "ArrowRight", preventDefault: function () {}});
  assert.ok(element("note-row").value.includes("[indices or NONE]"));
}
process.stdout.write("Offline viewer control contract passed\n");
